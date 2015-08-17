{-# LANGUAGE FlexibleInstances #-}
import Criterion.Main
import PCA.QModels
import PCA.Distribution
import Numeric.LinearAlgebra.Data (Vector, Matrix, fromRows, toRows)
import qualified Numeric.LinearAlgebra.HMatrix as H

import Control.DeepSeq

instance NFData (Distr Double Double) where
  rnf a = seq a ()

instance NFData (Distr (Vector Double) (Vector Double)) where
  rnf a = seq a ()

instance NFData (Distr (Vector Double) (Matrix Double)) where
  rnf a = seq a ()

instance NFData (Distr (Matrix Double) (Matrix Double)) where
  rnf a = seq a ()

prepareData n d = do
  trainT <- H.rand n d
  let tau = initQtau
      mu = initQmu d
      alpha = initQalpha d
      w = initQW d
      sigmaX = calcSigmaX tau w
  return (trainT, tau, mu, alpha, w, sigmaX)


foolCalcMX
    :: Matrix Double
    -> Distr Double Double
    -> Matrix Double
    -> Distr (Matrix Double) (Matrix Double)
    -> Distr (Vector Double) (Matrix Double)
    -> Matrix Double
foolCalcMX trainT tau sigmaX w mu =
    fromRows
        (map
             (\r ->
                   coef H.#> (r - mean mu))
             (toRows trainT))
  where
    coef = mean tau `H.scale` sigmaX H.<> H.tr (mean w)

main =
    defaultMain
        [ env (prepareData 100 3) $
          \small ->
               bgroup
                   "Small 100x3"
                   [ bench "Fool" $
                     whnf
                         (\(trainT,tau,mu,alpha,w,sigmaX) ->
                               foolCalcMX trainT tau sigmaX w mu)
                         small
                   , bench "Matrix" $
                     whnf
                         (\(trainT,tau,mu,alpha,w,sigmaX) ->
                               calcMX trainT tau sigmaX w mu)
                         small]
        , env (prepareData 10000 30) $
          \big ->
               bgroup
                   "Big 10000x30"
                   [ bench "Fool" $
                     whnf
                         (\(trainT,tau,mu,alpha,w,sigmaX) ->
                               foolCalcMX trainT tau sigmaX w mu)
                         big
                   , bench "Matrix" $
                     whnf
                         (\(trainT,tau,mu,alpha,w,sigmaX) ->
                               calcMX trainT tau sigmaX w mu)
                         big]]
