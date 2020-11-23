#ifndef __octree
#define __octree

#include <vector>
#include "create_octree.h"
using namespace std;


// 8分木のノードを表す構造体
// struct octreeNode {
//   vector <int> pInd;
//   double c[3];
//   octreeNode *cOctreeNode[2][2][2];
// };


// クラス octree の宣言
class octree {
private:
  vector <int> pInd;
public:
  octreeNode *octreeRoot;
  octree(float points[], size_t np, double range[], int nMin);
};

#endif
