import sys
import argparse
from wrapper import LinearRegression


def main():
  parser = argparse.ArgumentParser(description='Train linear regression')
  parser.add_argument("input_x", help=r"Path to file with feature values, delimiter='\n'", type=str)
  parser.add_argument("input_y", help=r"Path to file with target variable, delimiter='\n'", type=str)
  parser.add_argument("-i", "--iters", help="Number of gradient descent iterations", type=int, default=1000)
  parser.add_argument("-a", "--alpha", help="Learning rate", type=float, default=0.01)

  args = parser.parse_args()
  
  x_file_name = args.input_x
  y_file_name = args.input_y

  iters = args.iters
  alpha = args.alpha
  

  with open(x_file_name, 'r') as f:
    x = [float(line.strip()) for line in f]

  with open(y_file_name, 'r') as f:
    y = [float(line.strip()) for line in f]

  lr = LinearRegression()
  theta = lr.train(x, y, iters=iters, alpha=alpha)

  print(f'theta_0 = {theta[0]}')
  print(f'theta_1 = {theta[1]}')


if __name__ == '__main__':
  main()