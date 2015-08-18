
module PCA.QModels where

import PCA.Distribution
import Numeric.LinearAlgebra.Data
       (Matrix, Vector, ident, konst, size, toRows, fromRows, diag, diagl,
        takeDiag)
import qualified Numeric.LinearAlgebra.HMatrix as H


initQtau :: Distr Double Double
initQtau = Gamma 1.0e-3 1.0e-3

initQmu :: Int -> Distr (Vector Double) (Matrix Double)
initQmu d = MNormal (konst 0 d) (ident d)

initQalpha :: Int -> Distr (Vector Double) (Vector Double)
initQalpha d = MGamma (konst 1.0e-3 d) (konst 1.0e-3 d)

initQW :: Int -> Distr (Matrix Double) (Matrix Double)
initQW d = MatrixMNormal (ident d) (ident d)

initQX :: Matrix Double -> Distr (Matrix Double) (Matrix Double)
initQX trainT = let (_, d) = size trainT
                    sigmaX = calcSigmaX d initQtau (initQW d)
                    mX = calcMX trainT initQtau sigmaX (initQW d) (initQmu d)
                in MatrixMNormal mX sigmaX

type QModel = (Matrix Double, Vector Double)  -- ^ Mode of W and \mu

calculateQ :: Matrix Double -> QModel
calculateQ trainT =
    go 5 initQtau (initQmu d) (initQalpha d) (initQW d) (initQX trainT) trainT
  where
    (n,d) = size trainT
    go 0 _ mu _ w _ _ = (mean w, mean mu)
    go k (Gamma aTau bTau) mu (MGamma aAlpha bAlpha) w x t =
        go (k - 1) ntau nmu nalpha nw nx t
      where
        ntau = Gamma (calcAtau n d aTau) (calcBtau bTau t mu w x)
        nalpha = MGamma (calcAalpha d aAlpha) (calcBalpha bAlpha w)
        nSigmaMu = calcSigmaMu n ntau mu
        nmu = MNormal (calcMMu t ntau nSigmaMu w x) nSigmaMu
        nSigmaW = calcSigmaW ntau nalpha x
        nw = MatrixMNormal (calcMW ntau nSigmaW x t nmu) nSigmaW
        nSigmaX = calcSigmaX d ntau nw
        nx = MatrixMNormal (calcMX t ntau nSigmaX nw nmu) nSigmaX


calcSigmaX
    :: Int
    -> Distr Double Double
    -> Distr (Matrix Double) (Matrix Double)
    -> Matrix Double
calcSigmaX d tau w = H.inv (ident d + mean tau `H.scale` wSquareMean)
  where
    mW = mean w
    wSquareMean = variance w + (H.tr mW H.<> mW)

calcMX
    :: Matrix Double
    -> Distr Double Double
    -> Matrix Double
    -> Distr (Matrix Double) (Matrix Double)
    -> Distr (Vector Double) (Matrix Double)
    -> Matrix Double
calcMX trainT tau sigmaX w mu =
    H.tr
        (mean tau `H.scale` sigmaX H.<> H.tr (mean w) H.<>
         H.tr
             (fromRows
                  (map
                       (\r ->
                             r - mean mu)
                       (toRows trainT))))

calcSigmaMu
    :: Int
    -> Distr Double Double
    -> Distr (Vector Double) (Matrix Double)
    -> Matrix Double
calcSigmaMu n tau mu =
    vMu + 1 / (fromIntegral n * mean tau) `H.scale` ident d
  where
    vMu = variance mu
    (_,d) = size vMu

calcMMu
    :: Matrix Double
    -> Distr Double Double
    -> Matrix Double
    -> Distr (Matrix Double) (Matrix Double)
    -> Distr (Matrix Double) (Matrix Double)
    -> Vector Double
calcMMu trainT tau sigmaMu w x =
    (mean tau `H.scale` sigmaMu) H.#>
    foldr
         (flip (+))
         (konst 0 d)
         (toRows (trainT - H.tr (mW H.<> H.tr mX)))
  where
    mW = mean w
    mX = mean x
    (_, d) = size trainT


calcSigmaW
    :: Distr Double Double
    -> Distr (Vector Double) (Vector Double)
    -> Distr (Matrix Double) (Matrix Double)
    -> Matrix Double
calcSigmaW tau alpha x =
    H.inv
        (diag (mean alpha) +
         mean tau `H.scale`
         (fromIntegral n `H.scale` variance x + H.tr mX H.<> mX))
  where
    mX = mean x
    (n,_) = size mX

calcMW
    :: Distr Double Double
    -> Matrix Double
    -> Distr (Matrix Double) (Matrix Double)
    -> Matrix Double
    -> Distr (Vector Double) (Matrix Double)
    -> Matrix Double
calcMW tau sigmaW x trainT mu =
    (mean tau `H.scale` sigmaW) H.<>
    (H.tr mX H.<>
     (fromRows
          (map
               (\tn ->
                     tn - mMu)
               (toRows trainT))))
  where
    mX = mean x
    mMu = mean mu


calcAalpha :: Int -> Vector Double -> Vector Double
calcAalpha d = H.cmap (+ fromIntegral d / 2)

calcBalpha :: Vector Double
           -> Distr (Matrix Double) (Matrix Double)
           -> Vector Double
calcBalpha b w = b + 0.5 `H.scale` (takeDiag (H.tr mW H.<> mW))
  where
    mW = mean w

calcAtau :: Int -> Int -> Double -> Double
calcAtau n d a = a + fromIntegral (n * d) / 2


-- | TODO: Formula is not correct.
-- Skip part Tr (\langle W^TW \rangle \langle x_n x^T_n 'rangle')
-- because i think it is not correct in that result must be number
calcBtau
    :: Double
    -> Matrix Double
    -> Distr (Vector Double) (Matrix Double)
    -> Distr (Matrix Double) (Matrix Double)
    -> Distr (Matrix Double) (Matrix Double)
    -> Double
calcBtau b trainT mu w x =
    b +
    0.5 *
    (fromIntegral n * norm mMu +
     H.sumElements (takeDiag (trainT H.<> H.tr trainT)) +
     2 * (H.sumElements (H.tr (mW H.<> H.tr mX) H.#> mMu)) -
     2 * (H.sumElements (takeDiag (trainT H.<> mW H.<> H.tr mX))) -
     2 * (H.sumElements (trainT H.#> mMu)))
  where
    (n,_) = size trainT
    norm v = v `H.dot` v
    mMu = mean mu
    mW = mean w
    mX = mean x
