-- Cubic Bézier curve sampling for use with animation easing.
--
-- Code by Rafael Navega (2022).
-- License: Public Domain.
--
-- Main resource:
-- Numerical Recipes in C, "Chapter 5. Evaluation of Functions" (page 184)
--
-- References:
--
-- Synfig source code (function intersect_cubic):
--     https://github.com/synfig/synfig/blob/master/synfig-core/src/synfig/rendering/primitive/intersector.cpp#L311-L346
--
-- Blender source code (function findzero):
--     https://github.com/blender/blender/blob/master/source/blender/blenkernel/intern/fcurve.c#L1456
--
--     A control-point repositioning code that ensures the curve only has a single solution
--     vertically (function correct_bezpart):
--     https://github.com/martijnberger/blender/blob/master/source/blender/blenkernel/intern/fcurve.c#L1826-L1866
--
-- STB_Truetype source code (function stbtt__solve_cubic):
--     https://github.com/nothings/stb/blob/master/stb_truetype.h#L4543
--
-- FQS source code (function single_cubic):
--     https://github.com/NKrvavica/fqs/blob/master/fqs.py#L47
--
-- Articles:
-- https://trans4mind.com/personal_development/mathematics/polynomials/cubicAlgebra.htm
-- https://en.wikipedia.org/wiki/Cubic_equation#Depressed_cubic
-- https://pomax.github.io/bezierinfo/#extremities


local mathPow = math.pow
local mathSqrt = math.sqrt
local mathAcos = math.acos
local mathCos = math.cos


local EasingCurve = { }
EasingCurve.__index = EasingCurve


function EasingCurve.create(p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y)
    local c = {bezier=love.math.newBezierCurve(p0x, p0y, p1x, p1y, p2x, p2y, p3x, p3y)}
    setmetatable(c, EasingCurve)
    c:updateCache()
    return c
end


function EasingCurve:setControlPoint(index, x, y)
    self.bezier:setControlPoint(index, x, y)
    self:updateCache()
end


-- === function findTFromTime ===
-- NOTE: optimized version below, as findTFromTimeFast().
--
-- Parameters:
-- evalTime: the X coordinate of the vertical ray sampling the curve.
-- p0, p1, p2, p3: the X coordinates of the 4 control points of the cubic Bézier curve.
--
-- Returns:
-- The 't' parameter for the curve where the intersection happened, or returns
-- 'nil' in case there are multiple intersections (when the curve waves around the
-- vertical ray).
function EasingCurve:findTFromTime(evalTime, p0, p1, p2, p3)
    -- Take the "power basis" form of the Cubic Bézier:
    -- P(t) = (1 - t)^3 * P0 + 3t(1-t)^2 * P1 + 3t² (1-t) * P2 + t³ * P3

    -- We find the 4 coefficients from the full cubic polynomial. To do that,
    -- you run the above form through WolframAlpha to get the expanded version,
    -- then combine what you can to bring it down to this cubic polynomial form...
    -- (a t³ + b t² + c t + d = 0)
    --
    -- These a,b,c,d coefficients are based on the control points:
    -- (p3-p0+3*(p1-p2)) * t³ + (3*(p2-2*p1+p0)) * t² + 3*(p1-p0) * t + p0
    --
    -- To place the vertical ray at 'evalTime', subtract this value from
    -- the coefficient (d), moving the whole curve backwards so the ray is
    -- at the origin (x=0).
    local c3 = p3 + 3.0 * (p1 - p2) - p0
    local c2 = 3.0 * (p2 - 2.0 * p1 + p0)
    local c1 = 3.0 * (p1 - p0)
    local c0 = p0 - evalTime

    -- The root finding formula uses a reduced cubic polynomial, which is the
    -- full polynomial (ax³ + bx² + cx + d = 0) as divided by (a), giving:
    -- (x³ + ax² + bx + c = 0)
    --
    -- So these a,b,c coefficients below are from the *reduced* polynomial.
    local a = c2 / c3
    local b = c1 / c3
    local c = c0 / c3

    -- Neat optimization from the Blender source code: if you precalculate (a/3), the
    -- expressions below can be simplified.
    local aDiv3 = a / 3.0

    -- From the book, find values Q and R.

    -- Starting with the original formula for Q:
    -- Q = (a² - 3b) / 9
    --
    -- You expand and replace any (a/3) with the precomputed value:
    -- a²/9 - 3b/9 =
    -- a²/9 - b/3 =
    -- (a/3)² - b/3 =
    -- x² - b/3
    local Q = (aDiv3*aDiv3) - b/3.0

    -- Starting with the original formula for R:
    -- R = (2a³ - 9ab + 27c) / 54
    --
    -- The alternative/rearranged formula of which is (thanks WolframAlpha):
    -- R = a³/27 - ab/6 + c/2
    --
    -- You expand and replace any (a/3) as well:
    -- R = a³/27  - ab/6      + c/2 =
    --     (a/3)³ - (a/3)*b/2 + c/2 =
    --     x³     - xb/2      + c/2 =
    --     (2*x³   - x*b       + c) / 2
    local R = ((2 * aDiv3 * aDiv3 * aDiv3) - (aDiv3 * b) + c) / 2.0

    local RR = R * R
    local QQQ = Q * Q * Q
    if RR < QQQ then
        -- This happens in case there's three distinct real roots (the ray
        -- crosses the curve in 3 different points).
        do
            error('EasingCurve:findTFromTime() -> Multiple solutions for sampling the curve')
        end
        -- Anyway, the code to get the three roots (if they're within the [0, 1] range) is this.
        -- To avoid this route, run the control points of your curve(s) through the "control-point
        -- fixing code" mentioned in the references, at the start of this file.
        local sqrt_Q = mathSqrt(Q)
        local theta = mathAcos(R / mathSqrt(QQQ))
        local t1 = -2.0 * sqrt_Q * mathCos(theta / 3.0) - aDiv3
        local t2 = -2.0 * sqrt_Q * mathCos((theta + 2.0 * math.pi) / 3.0) - aDiv3
        local t3 = -2.0 * sqrt_Q * mathCos((theta - 2.0 * math.pi) / 3.0) - aDiv3
        return (t1 >= 0.0 and t1 <= 1.0) and t1 or nil,
               (t2 >= 0.0 and t2 <= 1.0) and t2 or nil,
               (t3 >= 0.0 and t3 <= 1.0) and t3 or nil
    else
        -- At least one real root, the other two are either complex numbers or the same
        -- real number.
        -- We only want the first real root (a single value) to be the result.
        local A = (
            (R > 0.0 and -mathPow(R + mathSqrt(RR-QQQ), 1.0/3.0)) or
             mathPow(-R + mathSqrt(RR-QQQ), 1.0/3.0)
        )
        local t1 = A + (A == 0.0 and 0.0 or Q/A) - aDiv3
        return ((t1 >= 0.0 and t1 <= 1.0) and t1 or nil), nil, nil
    end
end


-- Cache/precompute anything we can, to be used with findTFromTimeFast().
-- Each time the curve is modified this function needs to be called so these
-- cached values can be updated.
-- Most of the time this only needs to be called once after the curve data
-- is loaded, because easing curves are usually static in that they exist as
-- keyframes authored on software such as Blender or others.
function EasingCurve:updateCache()
    local p0x = self.bezier:getControlPoint(1)
    local p1x = self.bezier:getControlPoint(2)
    local p2x = self.bezier:getControlPoint(3)
    local p3x = self.bezier:getControlPoint(4)

    local c3 = p3x + 3.0 * (p1x - p2x) - p0x
    local c2 = 3.0 * (p2x - 2.0 * p1x + p0x)
    local c1 = 3.0 * (p1x - p0x)
    local a = c2 / c3
    local b = c1 / c3
    self.inv_c3 = (1.0 / c3) / 2.0

    local aDiv3 = a / 3.0
    self.aDiv3 = aDiv3
    -- Note how partialR includes a division by 2, included in inv_c3 above used with
    -- coefficient (c) inside findTFromTimeFast()).
    self.partialR = ((2.0 * aDiv3 * aDiv3 * aDiv3) - (aDiv3 * b)) / 2.0
    self.Q = (aDiv3 * aDiv3) - b / 3.0
    self.neg_2_sqrt_Q = -2.0 * mathSqrt(self.Q)
    self.QQQ = self.Q * self.Q * self.Q
    self.srqt_QQQ = mathSqrt(self.QQQ)
    self.p0 = p0x
end


function EasingCurve:findTFromTimeFast(evalTime)
    local c = (self.p0 - evalTime) * self.inv_c3
    local R = self.partialR + c
    local RR = R * R

    if RR < self.QQQ then
        do
            error('EasingCurve:findTFromTimeFast() -> Multiple solutions for sampling the curve')
        end
        local theta = mathAcos(R / self.srqt_QQQ)
        local t1 = self.neg_2_sqrt_Q * mathCos(theta / 3.0) - self.aDiv3
        local t2 = self.neg_2_sqrt_Q * mathCos((theta + 2.0 * 3.1415926535898) / 3.0) - self.aDiv3
        local t3 = self.neg_2_sqrt_Q * mathCos((theta - 2.0 * 3.1415926535898) / 3.0) - self.aDiv3
        return (t1 >= 0.0 and t1 <= 1.0) and t1 or nil,
               (t2 >= 0.0 and t2 <= 1.0) and t2 or nil,
               (t3 >= 0.0 and t3 <= 1.0) and t3 or nil
    else
        local A = (
            (R > 0.0 and -mathPow(R + mathSqrt(RR-self.QQQ), 1.0/3.0)) or
             mathPow(-R + mathSqrt(RR-self.QQQ), 1.0/3.0)
        )
        local t1 = A + (A == 0.0 and 0.0 or self.Q/A) - self.aDiv3
        return ((t1 >= 0.0 and t1 <= 1.0) and t1 or nil), nil, nil
    end
end

local curve = EasingCurve.create(1.0, 2.115, 202.708, 207.152, 254.712, -442.834, 501, -185.099)
curve.bezier:translate(100,500)
curve:updateCache()


function love.load()
    love.window.setTitle('Cubic Bézier Sampling')
end


function love.update(dt)
end


function love.draw()
    love.graphics.print('Move the mouse to sample the curve vertically. Press Esc to quit.', 10, 10)
    love.graphics.line(curve.bezier:render())

    local WIDTH = 5
    love.graphics.setColor(0.0, 0.7, 1.0)
    local x1, y1 = curve.bezier:getControlPoint(1)
    local x2, y2 = curve.bezier:getControlPoint(2)
    love.graphics.line(x1, y1, x2, y2)
    love.graphics.rectangle('line', x2-WIDTH, y2-WIDTH, WIDTH+WIDTH, WIDTH+WIDTH)
    x1, y1 = curve.bezier:getControlPoint(4)
    x2, y2 = curve.bezier:getControlPoint(3)
    love.graphics.line(x1, y1, x2, y2)
    love.graphics.rectangle('line', x2-WIDTH, y2-WIDTH, WIDTH+WIDTH, WIDTH+WIDTH)

    local p0x, p0y = curve.bezier:getControlPoint(1)
    local p1x, p1y = curve.bezier:getControlPoint(2)
    local p2x, p2y = curve.bezier:getControlPoint(3)
    local p3x, p3y = curve.bezier:getControlPoint(4)

    local mouseX = love.mouse.getX()
    local mouseY = love.mouse.getY()

    -- Draw a horizontal line at the mouse.
    --love.graphics.setColor(1.0, 1.0, 1.0, 0.3)
    --love.graphics.line(0, mouseY, love.graphics.getWidth(), mouseY)

    local tx1, tx2, tx3 = curve:findTFromTimeFast(mouseX)
    love.graphics.setColor(1.0, 0.0, 0.0, 1.0)
    love.graphics.line(mouseX, 0, mouseX, love.graphics.getHeight())
    love.graphics.setColor(1.0, 1.0, 0.0, 1.0)
    local best_tx = tx1 or tx2 or tx3
    if best_tx ~= nil then
        local x, y = curve.bezier:evaluate(best_tx)
        love.graphics.rectangle('line', x-WIDTH, y-WIDTH, WIDTH+WIDTH, WIDTH+WIDTH)
    end

    love.graphics.setColor(1.0, 1.0, 1.0, 1.0)
    love.graphics.print(tostring(tx1) .. ' ' .. tostring(tx2) .. ' ' .. tostring(tx3), 10, 30)
end


function love.keypressed(key)
    if key == 'escape' then
        love.event.quit()
    end
end
