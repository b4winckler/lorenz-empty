-- Solve paramaters for quadratic maps whose renormalization is empty.

module Main where

import Control.Monad (when)
import Data.List (foldl')
import System.Environment (getArgs, getProgName)
import System.Exit (exitFailure)
import System.IO (hPutStrLn, stderr)
import qualified Data.ByteString.Lazy as B
import qualified Data.Csv as Csv

type R = Double
type MapR = R -> R

type R2 = (R, R)
type MapR2 = R2 -> R2

type ReturnTimes = (Int, Int)


-- Quadratic Lorenz branches (q), its derivative wrt. x (dq), inverse (iq), and
-- derivatives wrt. u and v (duq, dvq)
q0, q1, dq0, dq1, iq0, iq1, duq0, duq1, dvq0, dvq1 :: R -> R -> R -> MapR
q0 c u _ x = u * (1 - ((c-x)/c)^2)
q1 c _ v x = 1 + v*(-1 + ((x-c)/(1-c))^2)
dq0 c u _ x = 2*u*(c-x)/c^2
dq1 c _ v x = 2*v*(x-c)/(1-c)^2
iq0 c u _ x = c * ( 1 - sqrt (1 - x/u) )
iq1 c _ v x = c + (1-c) * sqrt ((x-1+v)/v)
duq0 c _ _ x = 1 - (1-x/c)^2
duq1 c _ _ x = 0
dvq0 c _ _ x = 0
dvq1 c _ _ x = ((x-c)/(1-c))^2 - 1

-- Newton iteration in 1D, given f and its derivative df
newton1d :: MapR -> MapR -> R -> R
newton1d f df x0 = snd $ head $ dropWhile converging $ zip xs (tail xs)
  where
    xs = iterate step x0
    step x = x - f x / df x
    converging (x,y) = (x-y)^2 > 1e-12

-- Newton iteration in 2D, given f and its inverse Jacobian idf
newton2d :: MapR2 -> (R2 -> MapR2) -> R2 -> R2
newton2d f idf z0 = snd $ head $ dropWhile converging $ zip zs (tail zs)
  where
    zs = iterate step z0
    step z@(x,y) = let (dx,dy) = idf z (f z) in (x-dx, y-dy)
    converging ((x,y), (x',y')) = (x-x')^2 + (y-y')^2 > 1e-12

-- Given c and return times (a,b), find (u,v) such that quadratic Lorenz with
-- parameters (c,u,v) has (a,b)-renormalization which is empty.
solve_uv' :: ReturnTimes -> R -> R2
solve_uv' (a,b) c = newton2d f idf ((1+c)/2, 1-c/2)
  where
    f (u,v) = ( iterate (q1 c u v) u !! a - c
              , iterate (q0 c u v) (1-v) !! b - c )
    idf (u,v) (u',v') = ( (d11*u' - d01*v')/det
                        , (-d10*u' + d00*v')/det )
      where
        d00 = product $ map (dq1 c u v) uorb
        d01 = foldl' (\x y -> dq1 c u v y * x + dvq1 c u v y) 0 uorb
        d10 = foldl' (\x y -> dq0 c u v y * x + duq0 c u v y) 0 vorb
        d11 = negate $ product $ map (dq0 c u v) vorb
        det = d00*d11 - d01*d10
        uorb = take a $ iterate (q1 c u v) u
        vorb = take b $ iterate (q0 c u v) (1-v)

solve_p' b c u v = newton1d f df $ (c + iq0 c u v c)/2
  where
    f x = iterate (q1 c u v) (q0 c u v x) !! b - x
    df x = let xs = take b $ iterate (q1 c u v) (q0 c u v x)
           in product (map (dq1 c u v) xs) * dq0 c u v x - 1

solve_q' a c u v = newton1d f df $ (c + iq1 c u v c)/2
  where
    f x = iterate (q0 c u v) (q1 c u v x) !! a - x
    df x = let xs = take a $ iterate (q0 c u v) (q1 c u v x)
           in product (map (dq0 c u v) xs) * dq1 c u v x - 1

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

solve_all' (a,b) c = [1-v, iq0 c u v c, p, c, q, iq1 c u v c, u, dist0, dist1]
  where
    (u,v) = solve_uv' (a,b) c
    p = solve_p' b c u v
    q = solve_q' a c u v
    dist0 = log ( dq0 c u v (iq0 c u v p) / dq0 c u v (iq0 c u v q) )
    dist1 = log ( dq1 c u v (iq1 c u v q) / dq1 c u v (iq1 c u v p) )

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

logStrLn = hPutStrLn stderr

main = do
  as <- getArgs
  when (length as /= 3) $ do
    progName <- getProgName
    logStrLn $ "usage: " ++ progName ++ " A B N"
    logStrLn "    A       return time on right"
    logStrLn "    B       return time on left"
    logStrLn "    N       number of grid points"
    exitFailure

  let a = read $ as !! 0 :: Int
      b = read $ as !! 1 :: Int
      n = read $ as !! 2 :: Int
  putStrLn "vconj\txi1\tp\tc\tq\teta1\tu\tdist_left\tdist_right"
  B.putStr $ toCsv $ map (solve_all' (a,b)) $ grid n
