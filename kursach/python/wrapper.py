import os
import ctypes

lib_path = os.path.abspath('libtest.so')

class LinearRegression:
  def __init__(self):
    self.lr_dll = ctypes.CDLL(lib_path)

    self.lr_dll.train.restype = ctypes.POINTER(ctypes.c_double)
    self.lr_dll.train.argtypes = [
      ctypes.POINTER(ctypes.c_double),
      ctypes.POINTER(ctypes.c_double),
      ctypes.c_double,
      ctypes.c_int,
      ctypes.c_int
    ]

    self.lr_dll.predict_one.restype = ctypes.c_double
    self.lr_dll.predict_one.argtypes = [
      ctypes.c_double,
      ctypes.POINTER(ctypes.c_double)
    ]

    self.theta = None


  def train(self, x, y, iters=10, alpha=0.001):
    x_ = (ctypes.c_double * len(x))(*x)
    y_ = (ctypes.c_double * len(y))(*y)

    self.theta = self.lr_dll.train(x_, y_, alpha, iters, len(x))

    return self.theta[0], self.theta[1]

  def predict(self, x):
    if self.theta is None:
      raise Exception("Model isn't trained yet")

    return self.lr_dll.predict_one(x, self.theta)


if __name__ == '__main__':
  lr = LinearRegression()
    
  x = [1, 5, 10, 25, 33]
  y = [1, 4, 11, 23, 34]

  print("Train: theta = ", lr.train(x, y))
  
  print("Predict 15 = ", lr.predict(15))
