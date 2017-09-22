//
// Created by mmath on 7/10/17.
//
#include <cstdlib>
#include <iostream>

#include "RectangleScan.hpp"

int main(int argc, char* argv[]) {
  std::random_device rd;
  //std::uniform_int_distribution<int> dist(0, 1);
  //std::uniform_real_distribution<double> distXY(0, 1);
  std::vector<pyscan::Point<int>> points;

  size_t gSize = 10000;
  for (size_t i = 0; i < gSize; i++) {

      int red_blue = rand() % 2;
      float x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      float y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      points.push_back(pyscan::Point<int, 2>(red_blue, 1 - red_blue, x, y));
  }
  pyscan::Grid<int> grid(points.begin(), points.end(), 200);
  pyscan::Subgrid sg = pyscan::maxSubgridLinearG(grid, 200, .1, .1);
  std::cout << sg.toString() << std::endl;
  return 0;
}