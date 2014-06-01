-- Solve paramaters for quadratic maps whose renormalization is empty.

module Main where

import System.Environment (getArgs)
import qualified Data.Csv as Csv
import qualified Data.ByteString.Lazy as B

type EndoR = Double -> Double


-- Quadratic Lorenz branches (q), its derivative (dq) and inverse (iq)
q0, q1, dq0, dq1, iq0, iq1 :: Double -> Double -> Double -> EndoR
q0 c u v x = u * (1 - ((c-x)/c)^2)
q1 c u v x = 1 + v*(-1 + ((x-c)/(1-c))^2)
dq0 c u v x = 2*u*(c-x)/c^2
dq1 c u v x = 2*v*(x-c)/(1-c)^2
iq0 c u v x = c * ( 1 - sqrt (1 - x/u) )
iq1 c u v x = c + (1-c) * sqrt ((x-1+v)/v)

newton1d :: EndoR -> EndoR -> Double -> Double
newton1d f df x0 = snd $ head $ dropWhile converging $ zip xs (tail xs)
  where
    xs = iterate step x0
    step x = x - f x / df x
    converging (x,y) = (x-y)^2 > 1e-12

solve_p c u v = newton1d f df $ (c + iq0 c u v c)/2
  where
    f x = q1 c u v (q0 c u v x) - x
    df x = dq1 c u v (q0 c u v x) * dq0 c u v x - 1

solve_q c u v = newton1d f df $ (c + iq1 c u v c)/2
  where
    f x = q0 c u v (q1 c u v x) - x
    df x = dq0 c u v (q1 c u v x) * dq1 c u v x - 1

step_uv :: Double -> (Double, Double) -> (Double, Double)
step_uv c (u,v) = ( u - (d11*u' - d01*v')/det
                     , v - (-d10*u' + d00*v')/det )
  where
    d00 = 2*v*(u-c)
    d01 = (u-c)^2 - (1-c)^2
    d10 = c^2 - (c-1+v)^2
    d11 = -2*u*(c-1+v)
    det = d00*d11 - d01*d10
    u' = let c' = (1-c)^2 in c' + v * (-c' + (u - c)^2) - c*c'
    v' = let c' = c*c in u * (c' - (c - 1 + v)^2) - c*c'

solve_uv :: Double -> (Double, Double)
solve_uv c = snd $ head $ dropWhile converging $ zip xs (tail xs)
  where
    xs = iterate (step_uv c) ((1+c)/2, 1-c/2)
    converging ((x0,y0), (x1,y1)) = (x0-x1)^2 + (y0-y1)^2 > eps
    eps = 1e-12

grid :: Int -> [Double]
grid n = map ((/ fromIntegral (n+1)) . fromIntegral) [1..n]

dyadicGrid :: Int -> [Double]
dyadicGrid n = map ((2**) . negate . fromIntegral) [1..n]

solve_all c = [1-v, iq0 c u v c, p, c, q, iq1 c u v c, u, dist0, dist1]
  where
    (u,v) = solve_uv c
    p = solve_p c u v
    q = solve_q c u v
    dist0 = log ( dq0 c u v (iq0 c u v p) / dq0 c u v (iq0 c u v q) )
    dist1 = log ( dq1 c u v (iq1 c u v q) / dq1 c u v (iq1 c u v p) )

toCsv x = Csv.encodeWith encodeOpt x
  where
    encodeOpt = Csv.defaultEncodeOptions {
        Csv.encDelimiter = 9
      , Csv.encUseCrLf = False
      }

main = do
  as <- getArgs
  putStrLn "vconj\txi1\tp\tc\tq\teta1\tu\tdist_left\tdist_right"
  B.putStr $ toCsv $ map solve_all $ grid (read $ head as :: Int)
