
module QModels where

import Test.Tasty (TestTree, testGroup)
import Test.Tasty.HUnit
import PCA.QModels
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
                sigmaX = calcSigmaX tau w
            step "Test calculation of Q"
            assertEqual
                "Fool impl and matrix impl must give the same result"
                (calcMX trainT tau sigmaX w mu)
                (foolCalcMX trainT tau sigmaX w mu)
