{-# LANGUAGE FlexibleContexts #-}
module QModels where

import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit
import PCA.QModels
import PCA.Distribution
import Numeric.LinearAlgebra.Data (fromRows, toRows, size)
import qualified Numeric.LinearAlgebra.HMatrix as H


tests :: TestTree
tests = testGroup "QModels module tests" [unitTests]


unitTests :: TestTree
unitTests =
    testCaseSteps "Tests for Q model calculation" $
    \step ->
         do step "Preparing..."
            let n = 1000
                d = 30
            trainT <- H.rand n d
            let tau = initQtau
                mu = initQmu d
                alpha = initQalpha d
                w = initQW d
            step "Test calculation of Q"
            let sigmaX = calcSigmaX d tau w
            assertEqual "Size of Sigma_X must dxd" (size sigmaX) (d,d)
            let mX = calcMX trainT tau sigmaX w mu
            assertEqual
                "Fool impl and matrix impl must give the same result"
                mX
                (foolCalcMX trainT tau sigmaX w mu)
            let sigmaMu = calcSigmaMu n tau mu
            assertEqual "Size of Sigma_Mu is dxd" (size sigmaMu) (d,d)


foolCalcMX trainT tau sigmaX w mu =
    fromRows
        (map
             (\r ->
                   coef H.#> (r - mean mu))
             (toRows trainT))
  where
    coef = mean tau `H.scale` sigmaX H.<> H.tr (mean w)
