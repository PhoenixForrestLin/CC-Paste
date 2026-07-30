--=====================================================================
-- heli-fcs.lua | CC:Tweaked 直升机飞控（串级 PID）
-- 依赖：aero、sublevel API
-- 输出：主桨 / 左发 / 右发，均为 [-255, 255]
--=====================================================================

---------------------------------------------------------------------
-- 0. 配置
---------------------------------------------------------------------
local cfg = {
  loopHz = 20,                 -- 控制频率（CC 一 tick = 0.05s，20Hz 刚好）

  -- 期望值上限（外环输出的限幅）
  maxClimbRate  = 6.0,         -- m/s
  maxCruiseSpd  = 12.0,        -- m/s
  maxYawRate    = 1.2,         -- rad/s

  -- 摇杆 → 期望偏移速率（松杆即停，setpoint 不再移动）
  stick = {
    shiftRate = 8.0,           -- 前后杆满舵时目标点每秒前移 m
    yawRate   = 1.0,           -- 转向杆满舵时目标航向每秒转 rad
    climbRate = 5.0,           -- 上下档位 15 档对应的目标高度变化 m/s
    leash     = 60.0,          -- 摇杆最多把目标点拉离原任务点多少 m
    deadband  = 0.02,
  },

  -- 主桨推力标定：单位输出产生的推力（N/unit），先粗标，之后由 hoverTrim 自学习
  -- 实测 mass = 115, g = 11 → 悬停总需推力 ≈ 1265 N
  thrustPerUnitMain = 10.54,     -- N / unit
  hoverTrimTau      = 4.0,     -- 悬停配平学习时间常数 s

  -- 环境兜底值（实测 aero.getDefault(): gravity = 0,-11,0 / pressure = 1 / universalDrag = 0.09）
  gravityFallback   = 11.0,
  pressureFallback  = 1.0,
  massFallback      = 115.0,

  -- 到达判据
  arriveRadius = 2.0,          -- m
  arriveAlt    = 0.7,          -- m
  headingLock  = 0.35,         -- rad，机头误差小于此值才开始前推
  faceDeadzone = 3.0,          -- m，距离小于此值不再强行对准机头

  -- 符号约定（若实机反向，把对应项改为 -1）
  sign = { yaw = -1, fwd = -1, main = -1 },

  -- 输出限制
  slewMain = 900,              -- 每秒最大变化量（unit/s），防止硬阶跃
  slewEng  = 1200,

  -- 失效保护
  maxTiltDot   = 0.55,         -- 机体上向量与世界上向量点积低于此值判定失姿
  poseTimeout  = 0.6,          -- s，读不到位姿即进入 FAILSAFE
}

-- 增益：{kp, ki, kd}
local gains = {
  alt      = { kp = 0.90, ki = 0.00, kd = 0.15 },  -- 高度 → 期望垂速
  vz       = { kp = 28.0, ki = 12.0, kd = 3.0  },  -- 垂速 → 主桨增量
  pos      = { kp = 0.55, ki = 0.00, kd = 0.10 },  -- 水平距离 → 期望前速
  vx       = { kp = 22.0, ki = 8.00, kd = 2.0  },  -- 前速 → 共通油门
  heading  = { kp = 1.60, ki = 0.00, kd = 0.05 },  -- 航向 → 期望偏航角速度
  yawRate  = { kp = 60.0, ki = 20.0, kd = 4.0  },  -- 偏航角速度 → 差量
}

--设备配置

local dashPanel = peripheral.wrap("control_panel_2")

local device = {
  -- 主桨
  mainRotor = peripheral.wrap("Create_RotationSpeedController_12"),

  -- 左右发
  leftEngine = peripheral.wrap("Create_RotationSpeedController_10"),
  rightEngine = peripheral.wrap("Create_RotationSpeedController_9"),

  -- 输入设备
  switch = dashPanel.getModule("switch"),
  controlLever = dashPanel.getModule("control_lever"),
  joystick = dashPanel.getModule("joystick"),
}

---------------------------------------------------------------------
-- 1. 数学工具
---------------------------------------------------------------------
local function clamp(v, lo, hi) if v < lo then return lo elseif v > hi then return hi end return v end
local function dz(v, d) if math.abs(v) < d then return 0 end return v end
local function wrapPi(a)
  while a >  math.pi do a = a - 2 * math.pi end
  while a < -math.pi do a = a + 2 * math.pi end
  return a
end

-- 把任意形状（{x,y,z} 或 {[1],[2],[3]} 或 vector）统一成 {x,y,z}
local function toVec(t)
  if type(t) ~= "table" then return nil end
  if t.x ~= nil then return { x = t.x, y = t.y or 0, z = t.z or 0 } end
  return { x = t[1] or 0, y = t[2] or 0, z = t[3] or 0 }
end

-- 四元数旋转向量：v' = q * v * q^-1
local function qRot(q, v)
  local x, y, z, w = q.x, q.y, q.z, q.w
  local tx = 2 * (y * v.z - z * v.y)
  local ty = 2 * (z * v.x - x * v.z)
  local tz = 2 * (x * v.y - y * v.x)
  return {
    x = v.x + w * tx + (y * tz - z * ty),
    y = v.y + w * ty + (z * tx - x * tz),
    z = v.z + w * tz + (x * ty - y * tx),
  }
end

---------------------------------------------------------------------
-- 2. PID（测量微分 + 一阶滤波 + 反算抗饱和）
---------------------------------------------------------------------
local PID = {}
PID.__index = PID

function PID.new(g, outMin, outMax, dTau)
  return setmetatable({
    kp = g.kp, ki = g.ki, kd = g.kd,
    outMin = outMin, outMax = outMax,
    dTau = dTau or 0.08,
    i = 0, dFilt = 0, prev = nil,
  }, PID)
end

function PID:reset(keepI)
  if not keepI then self.i = 0 end
  self.dFilt, self.prev = 0, nil
end

-- ff: 前馈项（如悬停配平），直接加进输出并参与抗饱和
function PID:step(sp, meas, dt, ff)
  local err = sp - meas
  local p = self.kp * err

  local d = 0
  if self.prev then
    local raw = -(meas - self.prev) / dt          -- 对测量微分，避免设定值突跳
    local a = dt / (self.dTau + dt)
    self.dFilt = self.dFilt + a * (raw - self.dFilt)
    d = self.kd * self.dFilt
  end
  self.prev = meas

  self.i = self.i + self.ki * err * dt
  local raw = p + d + (ff or 0) + self.i
  local out = clamp(raw, self.outMin, self.outMax)
  if out ~= raw then self.i = self.i + (out - raw) end  -- 反算：积分不再继续膨胀
  return out, err
end

---------------------------------------------------------------------
-- 3. 传感器抽象（aero / sublevel）
---------------------------------------------------------------------
local sensor = {}

---------------------------------------------------------------------
-- 实测接口形状（控制台验证）：
--   sublevel.getLogicalPose()    → { orientation = w+xi+yj+zk, position = vec,
--                                    rotationPoint = vec, scale = vec }
--   sublevel.getLinearVelocity() → 多返回值 x, y, z（不是 table！）
--   sublevel.getAngularVelocity()→ 多返回值 x, y, z
--   sublevel.getMass()           → number（本机 115）
--   aero.getDefault()            → table：gravity = 0,-11,0 / pressure = 1 /
--                                  universalDrag = 0.09（不是 number！）
---------------------------------------------------------------------

local env = nil
local function getEnv()
  if env then return env end
  env = { gravity = cfg.gravityFallback, pressure = cfg.pressureFallback, drag = 0 }
  local ok, d = pcall(aero.getDefault)
  if ok and type(d) == "table" then
    local g = toVec(d.gravity)
    if g then
      local mag = math.sqrt(g.x * g.x + g.y * g.y + g.z * g.z)
      if mag > 1e-3 then env.gravity = mag end
    end
    env.pressure = tonumber(d.pressure) or env.pressure
    env.drag     = tonumber(d.universalDrag) or 0
  end
  return env
end

-- 多返回值 (x, y, z) → {x,y,z}；同时兼容返回 vector / table 的实现
local function readVec3(fn)
  if type(fn) ~= "function" then return nil end
  local ok, a, b, c = pcall(fn)
  if not ok then return nil end
  if type(a) == "number" then return { x = a, y = b or 0, z = c or 0 } end
  return toVec(a)
end

-- 四元数提取：优先 .w/.x/.y/.z；数组形式按打印顺序 w, x, y, z 解析
local function toQuat(rot)
  local q = { x = 0, y = 0, z = 0, w = 1 }
  if rot == nil then return q end
  local okw, w = pcall(function() return rot.w end)
  if okw and type(w) == "number" then
    q.w = w
    q.x = tonumber(rot.x) or 0
    q.y = tonumber(rot.y) or 0
    q.z = tonumber(rot.z) or 0
    return q
  end
  local ok1, first = pcall(function() return rot[1] end)
  if ok1 and type(first) == "number" then
    q.w = first
    q.x = tonumber(rot[2]) or 0
    q.y = tonumber(rot[3]) or 0
    q.z = tonumber(rot[4]) or 0
  end
  return q
end

local function readPose()
  local ok, pose = pcall(sublevel.getLogicalPose)
  if not ok or type(pose) ~= "table" then
    if type(sublevel.getLastPose) ~= "function" then return nil end
    ok, pose = pcall(sublevel.getLastPose)
    if not ok or type(pose) ~= "table" then return nil end
  end
  local pos = toVec(pose.position)
  if not pos then return nil end
  return { pos = pos, q = toQuat(pose.orientation or pose.rotation) }
end

-- 读取一帧完整状态
function sensor.read()
  local pose = readPose()
  if not pose then return nil end

  local fwd = qRot(pose.q, { x = 0, y = 0, z = 1 })   -- 机体前向（若实机为 -Z，改这里）
  local up  = qRot(pose.q, { x = 0, y = 1, z = 0 })

  local hLen = math.sqrt(fwd.x * fwd.x + fwd.z * fwd.z)
  if hLen < 1e-6 then hLen = 1e-6 end
  local fh = { x = fwd.x / hLen, z = fwd.z / hLen }    -- 水平前向单位向量

  local vel = readVec3(sublevel.getLinearVelocity)  or { x = 0, y = 0, z = 0 }
  local omg = readVec3(sublevel.getAngularVelocity) or { x = 0, y = 0, z = 0 }
  local e = getEnv()

  return {
    pos     = pose.pos,
    up      = up,
    fh      = fh,
    heading = math.atan2(fh.x, fh.z),
    vFwd    = (vel.x * fh.x + vel.z * fh.z),            -- 机体前向速度
    vLat    = (vel.x * fh.z - vel.z * fh.x),            -- 机体侧向速度（仅监视）
    vz      = vel.y,
    yawRate = cfg.sign.yaw * dz(omg.y, 1e-4),           -- 静止读数约 1e-20，直接当 0
    mass    = tonumber(sublevel.getMass and sublevel.getMass()) or cfg.massFallback,
    gravity = e.gravity,
    drag    = e.drag,
  }
end

-- 空气密度补偿：基准气压取 aero.getDefault().pressure（实测 = 1），不是 getDefault() 本身
local function airFactor(pos)
  if type(aero.getAirPressure) ~= "function" then return 1 end
  local ok, p = pcall(aero.getAirPressure, vector.new(pos.x, pos.y, pos.z))
  if not ok or type(p) ~= "number" then
    ok, p = pcall(aero.getAirPressure, pos.x, pos.y, pos.z)   -- 兼容多参数签名
  end
  if not ok or type(p) ~= "number" or p <= 0 then return 1 end
  local base = getEnv().pressure
  return clamp(base / p, 0.6, 2.0)
end

---------------------------------------------------------------------
-- 4. 输入抽象（替换成你的实际来源：rednet / 触摸屏 / 红石接口）
---------------------------------------------------------------------
local stickState = { 
  pitch = 0, 
  yaw = 0, 
  lift = device.controlLever.getValue(), 
  liftUp =  device.switch.getState()
}

local function readSticks()
  -- TODO: 接入实际输入源，只需保证：
  --   pitch : [-1,1] 精度 0.01，正 = 向前
  --   yaw   : [-1,1] 精度 0.01，正 = 右转
  --   lift  : 整数 [0,15]
  --   liftUp: bool，true = 上升
  local s = stickState
  return {
    pitch  = clamp(dz(s.pitch or 0, cfg.stick.deadband), -1, 1),
    yaw    = clamp(dz(s.yaw   or 0, cfg.stick.deadband), -1, 1),
    lift   = clamp(math.floor(s.lift or 0), 0, 15) / 15,
    liftUp = s.liftUp ~= false,
  }
end

-- 后台监听输入（示例：rednet 广播 {pitch=..,yaw=..,lift=..,up=..}）
local function inputTask()
  while true do
    local _, msg = rednet.receive("heli_stick")
    if type(msg) == "table" then
      stickState.pitch  = tonumber(msg.pitch) or 0
      stickState.yaw    = tonumber(msg.yaw) or 0
      stickState.lift   = tonumber(msg.lift) or 0
      stickState.liftUp = msg.up ~= false
    end
  end
end

---------------------------------------------------------------------
-- 5. 输出（替换成你的实际执行器接口）
---------------------------------------------------------------------
local outState = {
  main = device.mainRotor.getTargetSpeed(),
  left = device.leftEngine.getTargetSpeed(),
  right = device.rightEngine.getTargetSpeed()
}

local function writeOutputs(main, left, right)
  outState.main, outState.left, outState.right = main, left, right
  device.mainRotor.setTargetSpeed(outState.main)
  device.leftEngine.setTargetSpeed(outState.left)
  device.rightEngine.setTargetSpeed(outState.right)
end

---------------------------------------------------------------------
-- 6. 飞控主体
---------------------------------------------------------------------
local FCS = {}
FCS.__index = FCS

function FCS.new(tx, ty, tz, cruiseAlt)
  local self = setmetatable({}, FCS)
  self.mission = { x = tx, y = ty, z = tz, cruise = cruiseAlt }
  self.sp = { x = nil, z = nil, alt = cruiseAlt, hdg = nil }  -- 实际执行的设定值
  self.offset = { x = 0, z = 0, alt = 0, hdg = 0 }            -- 摇杆偏移量
  self.phase = "INIT"
  self.hoverTrim = nil
  self.out = { main = 0, left = 0, right = 0 }
  self.lostFor = 0

  self.pidAlt = PID.new(gains.alt,     -cfg.maxClimbRate, cfg.maxClimbRate)
  self.pidVz  = PID.new(gains.vz,      -255, 255)
  self.pidPos = PID.new(gains.pos,     -cfg.maxCruiseSpd * 0.35, cfg.maxCruiseSpd)
  self.pidVx  = PID.new(gains.vx,      -255, 255)
  self.pidHdg = PID.new(gains.heading, -cfg.maxYawRate, cfg.maxYawRate)
  self.pidYaw = PID.new(gains.yawRate, -255, 255)
  return self
end

-- 摇杆只移动设定值：松杆 → 设定值冻结 → 原地悬停
function FCS:applySticks(st, s, dt)
  local o = self.offset
  o.x = o.x + cfg.stick.shiftRate * st.pitch * cfg.sign.fwd * s.fh.x * dt
  o.z = o.z + cfg.stick.shiftRate * st.pitch * cfg.sign.fwd * s.fh.z * dt
  o.hdg = wrapPi(o.hdg + cfg.stick.yawRate * st.yaw * dt)
  if st.lift > 0 then
    o.alt = o.alt + cfg.stick.climbRate * st.lift * (st.liftUp and 1 or -1) * dt
  end

  local r = math.sqrt(o.x * o.x + o.z * o.z)
  if r > cfg.stick.leash then
    o.x, o.z = o.x * cfg.stick.leash / r, o.z * cfg.stick.leash / r
  end
  o.alt = clamp(o.alt, -cfg.stick.leash, cfg.stick.leash)
end

function FCS:updatePhase(s)
  local m = self.mission
  local dx, dz_ = m.x - s.pos.x, m.z - s.pos.z
  local dist = math.sqrt(dx * dx + dz_ * dz_)

  if self.phase == "INIT" then
    self.phase = (math.abs(s.pos.y - m.cruise) > cfg.arriveAlt) and "TAKEOFF" or "CRUISE"
  elseif self.phase == "TAKEOFF" then
    if math.abs(s.pos.y - m.cruise) < cfg.arriveAlt then self.phase = "CRUISE" end
  elseif self.phase == "CRUISE" then
    if dist < cfg.arriveRadius then self.phase = "APPROACH" end
  elseif self.phase == "APPROACH" then
    if dist < cfg.arriveRadius and math.abs(s.pos.y - m.y) < cfg.arriveAlt then
      self.phase = "HOLD"
    elseif dist > cfg.arriveRadius * 2.5 then
      self.phase = "CRUISE"
    end
  end
  return dist
end

function FCS:step(s, st, dt)
  local m = self.mission
  self:applySticks(st, s, dt)
  local dist = self:updatePhase(s)

  -- === 设定值装配 ===
  local tgtX = m.x + self.offset.x
  local tgtZ = m.z + self.offset.z
  local baseAlt = (self.phase == "TAKEOFF" or self.phase == "CRUISE") and m.cruise or m.y
  local tgtAlt = baseAlt + self.offset.alt

  local ex, ez = tgtX - s.pos.x, tgtZ - s.pos.z
  local hDist = math.sqrt(ex * ex + ez * ez)
  local bearing = (hDist > 1e-3) and math.atan2(ex, ez) or s.heading
  local tgtHdg = wrapPi(((hDist > cfg.faceDeadzone) and bearing or s.heading) + self.offset.hdg)

  -- === 垂直通道（串级）===
  local vzCmd = self.pidAlt:step(tgtAlt, s.pos.y, dt)
  if st.lift > 0 then                              -- 有升降杆时直接给速度指令，手感更硬
    vzCmd = cfg.maxClimbRate * st.lift * (st.liftUp and 1 or -1)
  end
  -- 悬停前馈：需要推力 = m·g（实测 115 × 11 ≈ 1265 N），除以单位输出推力得到 unit 数
  local hoverFF = (s.mass * s.gravity)
                  / math.max(cfg.thrustPerUnitMain, 1e-6) * airFactor(s.pos)
  hoverFF = clamp(self.hoverTrim or hoverFF, 0, 255)
  local main = self.pidVz:step(vzCmd, s.vz, dt, hoverFF)

  -- 悬停配平自学习：稳定悬停时把当前主桨值滤成新的前馈基准
  if math.abs(s.vz) < 0.25 and math.abs(tgtAlt - s.pos.y) < 0.8 then
    local a = dt / (cfg.hoverTrimTau + dt)
    self.hoverTrim = (self.hoverTrim or main) + a * (main - (self.hoverTrim or main))
  end

  -- === 航向通道（串级）===
  local hdgErr = wrapPi(tgtHdg - s.heading)
  local yawRateCmd = self.pidHdg:step(0, -hdgErr, dt)
  local diff = self.pidYaw:step(yawRateCmd, s.yawRate, dt)

  -- === 水平通道（串级）===
  local vFwdCmd = 0
  if hDist > 0.3 then
    local along = self.pidPos:step(hDist, 0, dt)   -- 距离 → 期望前速
    -- 机头没对准时按 cos 衰减，超过阈值几乎不前推，先转再走
    local gate = math.cos(clamp(math.abs(hdgErr), 0, math.pi / 2))
    if math.abs(hdgErr) > cfg.headingLock and hDist > cfg.faceDeadzone then gate = gate * 0.15 end
    vFwdCmd = clamp(along * gate, -cfg.maxCruiseSpd * 0.35, cfg.maxCruiseSpd)
  else
    self.pidPos:reset(true)
  end
  local common = self.pidVx:step(vFwdCmd, s.vFwd, dt) * cfg.sign.fwd

  -- === 混控 ===
  local left, right = common - diff, common + diff
  local peak = math.max(math.abs(left), math.abs(right))
  if peak > 255 then                                -- 等比缩放，保住差速比例（转向优先级不丢）
    left, right = left * 255 / peak, right * 255 / peak
  end

  -- === 限幅 + 变化率限制 ===
  local function slew(prev, want, rate)
    local d = clamp(want - prev, -rate * dt, rate * dt)
    return clamp(prev + d, -255, 255)
  end
  self.out.main  = slew(self.out.main,  main * cfg.sign.main, cfg.slewMain)
  self.out.left  = slew(self.out.left,  left,  cfg.slewEng)
  self.out.right = slew(self.out.right, right, cfg.slewEng)

  writeOutputs(
    math.floor(self.out.main + 0.5),
    math.floor(self.out.left + 0.5),
    math.floor(self.out.right + 0.5)
  )

  self.debug = {
    phase = self.phase, dist = hDist, altErr = tgtAlt - s.pos.y,
    hdgErr = hdgErr, vzCmd = vzCmd, vFwdCmd = vFwdCmd, trim = self.hoverTrim,
  }
end

-- 失姿/失读保护：主桨维持配平缓降，水平推力归零
function FCS:failsafe(dt)
  local main = clamp((self.hoverTrim or 0) * 0.85, 0, 255)
  self.out.main = self.out.main + clamp(main - self.out.main, -cfg.slewMain * dt, cfg.slewMain * dt)
  self.out.left, self.out.right = 0, 0
  self.pidVz:reset(true); self.pidVx:reset(); self.pidYaw:reset(); self.pidPos:reset(); self.pidHdg:reset()
  writeOutputs(math.floor(self.out.main + 0.5), 0, 0)
end

---------------------------------------------------------------------
-- 7. 主循环
---------------------------------------------------------------------
local function controlTask(fcs)
  local dtNom = 1 / cfg.loopHz
  local last = os.epoch("utc") / 1000
  fcs.ticks = 0
  while true do
    sleep(dtNom)                       -- 不再手动管理 timer id
    local now = os.epoch("utc") / 1000
    local dt = clamp(now - last, 0.01, 0.25)
    last = now
    fcs.ticks = fcs.ticks + 1

    local ok, err = pcall(function()
      local s = sensor.read()
      local bad = (s == nil) or (s.up.y < cfg.maxTiltDot)
      if bad then
        fcs.lostFor = fcs.lostFor + dt
        if fcs.lostFor > cfg.poseTimeout then fcs.phase = "FAILSAFE" end
        fcs:failsafe(dt)
      else
        fcs.lostFor = 0
        if fcs.phase == "FAILSAFE" then fcs.phase = "INIT" end
        fcs:step(s, readSticks(), dt)
      end
    end)
    if not ok then fcs.lastErr = tostring(err) end
  end
end

local function hudTask(fcs)
  while true do
    local d = fcs.debug
    if d then
      term.clear(); term.setCursorPos(1, 1)
      print(("PHASE %s"):format(d.phase))
      print(("dist %.2f  altErr %+.2f  hdgErr %+.2f"):format(d.dist, d.altErr, d.hdgErr))
      print(("vzCmd %+.2f  vFwdCmd %+.2f  trim %.1f"):format(d.vzCmd, d.vFwdCmd, d.trim or 0))
      print(("OUT main %d  L %d  R %d"):format(outState.main, outState.left, outState.right))
      print(("tick %d  dt %s"):format(fcs.ticks or 0, fcs.lastErr and "ERR" or "ok"))
      if fcs.lastErr then print(("ERR %s"):format(fcs.lastErr)) end
    end
    sleep(0.2)
  end
end

---------------------------------------------------------------------
-- 入口：flyTo(x, y, z, cruiseAlt)
---------------------------------------------------------------------
local function flyTo(x, y, z, cruiseAlt)
  local fcs = FCS.new(x, y, z, cruiseAlt or (y + 20))
  local tasks = { function() controlTask(fcs) end, function() hudTask(fcs) end }
  if rednet and peripheral.find("modem") then
    peripheral.find("modem", function(n) rednet.open(n) end)
    table.insert(tasks, inputTask)
  end
  parallel.waitForAny(table.unpack(tasks))
end

local a = { ... }
if #a >= 3 then
  flyTo(tonumber(a[1]), tonumber(a[2]), tonumber(a[3]), tonumber(a[4]))
else
  print("usage: heli-fcs <x> <y> <z> [cruiseAlt]")
end

return { flyTo = flyTo, FCS = FCS, cfg = cfg, gains = gains }
