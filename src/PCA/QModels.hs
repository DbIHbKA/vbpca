
module PCA.QModels where

import PCA.Distribution
import Numeric.LinearAlgebra.Data
       (Matrix, Vector, ident, konst, size, toRows, fromRows)
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
                    sigmaX = calcSigmaX initQtau (initQW d)
                    mX = calcMX trainT initQtau sigmaX (initQW d) (initQmu d)
                in MatrixMNormal mX sigmaX

type QModel = (Matrix Double, Vector Double)  -- ^ Mode of W and \mu

calculateQ :: Matrix Double -> QModel
calculateQ trainT = undefined
    -- go 20 initQtau (initQmu d) (initQalpha d) (initQW d) (initQX trainT)


calcSigmaX :: Distr Double Double
           -> Distr (Matrix Double) (Matrix Double)
           -> Matrix Double
calcSigmaX tau w = H.inv (ident d + mean tau `H.scale` wSquareMean)
  where
    mW = mean w
    wSquareMean = variance w + (H.tr mW H.<> mW)
    (_,d) = size mW

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

