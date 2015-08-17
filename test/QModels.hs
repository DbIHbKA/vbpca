{-# LANGUAGE FlexibleContexts #-}
module QModels where

import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit
import PCA.QModels
import PCA.Distribution
import Numeric.LinearAlgebra.Data
       (fromRows, toRows, size, (!), konst)
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
            assertEqual "Size of Sigma_X must dxd" (size sigmaX) (d, d)
            let mX = calcMX trainT tau sigmaX w mu
            assertEqual
                "Fool impl and matrix impl must give the same result"
                mX
                (foolCalcMX trainT tau sigmaX w mu)
            assertEqual "Size of mean X is Nxd" (size mX) (n,d)
            let x = MatrixMNormal mX sigmaX
            let sigmaMu = calcSigmaMu n tau mu
            assertEqual "Size of Sigma_Mu is dxd" (size sigmaMu) (d, d)
            let mMu = calcMMu trainT tau sigmaMu w x
            assertEqual "Size of m_Mu is 1xd" (size mMu) d
            assertEqual
                "Fool impl and matrix impl must give the same result"
                mMu
                (foolCalcMMu trainT tau sigmaMu w x)


foolCalcMX trainT tau sigmaX w mu =
    fromRows
        (map
             (\r ->
                   coef H.#> (r - mean mu))
             (toRows trainT))
  where
    coef = mean tau `H.scale` sigmaX H.<> H.tr (mean w)

foolCalcMMu trainT tau sigmaMu w x =
    (mean tau `H.scale` sigmaMu) H.#>
    foldr
         (\r resV ->
               resV + ((trainT ! r) - mW H.#> (mX ! r)))
         (konst 0 d)
         [0 .. (n - 1)]
  where
    mW = mean w
    mX = mean x
    (n,d) = size trainT
