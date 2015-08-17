{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE TypeFamilies #-}
module PCA.Distribution where

import Numeric.LinearAlgebra.Data (Vector, Matrix)


data family Distr m v

data instance Distr Double Double = Gamma Double Double

data instance
     Distr (Vector Double) (Vector Double) = MGamma (Vector Double)
                                                    (Vector Double)


data instance
     Distr (Vector Double) (Matrix Double) = MNormal (Vector Double)
                                                     (Matrix Double)


data instance
     Distr (Matrix Double) (Matrix Double) = MatrixMNormal (Matrix
                                                              Double)
                                                           (Matrix Double)

class Mean m v where
    mean :: Distr m v -> m

instance Mean Double Double where
    mean (Gamma s r) = s / r

instance Mean (Vector Double) (Vector Double) where
    mean (MGamma s r) = s / r

instance Mean (Vector Double) (Matrix Double) where
    mean (MNormal m v) = m

instance Mean (Matrix Double) (Matrix Double) where
    mean (MatrixMNormal m v) = m


class Variance m v  where
    variance :: Distr m v -> v

instance Variance (Vector Double) (Matrix Double) where
    variance (MNormal m v) = v
instance Variance (Matrix Double) (Matrix Double) where
    variance (MatrixMNormal m v) = v
