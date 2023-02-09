// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>

#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/InitialElementIds.hpp"
#include "Domain/Structure/SegmentId.hpp"
#include "Domain/ZCurve.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"

namespace {
// Test the Z-curve index computed for ElementIds in 1D. Tests
// domain::z_curve_index_from_element_id.
void test_z_curve_index_from_element_id_1d() {
  // The Z-curve does not depend on the block ID or grid index, so the choice of
  // these two values for this test is arbitrary
  const size_t block_id = 0;
  const size_t grid_index = 0;

  // Create SegmentIds from refinement 0 to 2
  // l0i0 refers to refinement 0, at index 0
  SegmentId l0i0(0, 0);

  SegmentId l1i0(1, 0);
  SegmentId l1i1(1, 1);

  SegmentId l2i0(2, 0);
  SegmentId l2i1(2, 1);
  SegmentId l2i2(2, 2);
  SegmentId l2i3(2, 3);

  // Create ElementIds with 1D SegmentIds from refinement 0 to 2
  const ElementId<1> e_l0x0(block_id, make_array<1>(l0i0), grid_index);

  const ElementId<1> e_l1x0(block_id, make_array<1>(l1i0), grid_index);
  const ElementId<1> e_l1x1(block_id, make_array<1>(l1i1), grid_index);

  const ElementId<1> e_l2x0(block_id, make_array<1>(l2i0), grid_index);
  const ElementId<1> e_l2x1(block_id, make_array<1>(l2i1), grid_index);
  const ElementId<1> e_l2x2(block_id, make_array<1>(l2i2), grid_index);
  const ElementId<1> e_l2x3(block_id, make_array<1>(l2i3), grid_index);

  // Check that ElementIds are mapped to their correct Z-curve index
  CHECK(domain::z_curve_index_from_element_id(e_l0x0) == 0);

  CHECK(domain::z_curve_index_from_element_id(e_l1x0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1) == 1);

  CHECK(domain::z_curve_index_from_element_id(e_l2x0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3) == 3);
}

// Test the Z-curve index computed for ElementIds in 2D. Tests
// domain::z_curve_index_from_element_id.
void test_z_curve_index_from_element_id_2d() {
  // The Z-curve does not depend on the block ID or grid index, so the choice of
  // these two values for this test is arbitrary
  const size_t block_id = 0;
  const size_t grid_index = 0;

  // Create SegmentIds from refinement 0 to 2
  // l0i0 refers to refinement 0, at index 0
  SegmentId l0i0(0, 0);

  SegmentId l1i0(1, 0);
  SegmentId l1i1(1, 1);

  SegmentId l2i0(2, 0);
  SegmentId l2i1(2, 1);
  SegmentId l2i2(2, 2);
  SegmentId l2i3(2, 3);

  // Create ElementIds with 2D SegmentIds from refinement 0 to 2 and check that
  // ElementIds are mapped to their correct Z-curve index
  const ElementId<2> e_l0x0_l0y0(block_id, make_array(l0i0, l0i0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l0y0) == 0);

  const ElementId<2> e_l1x0_l0y0(block_id, make_array(l1i0, l0i0), grid_index);
  const ElementId<2> e_l1x1_l0y0(block_id, make_array(l1i1, l0i0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l0y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l0y0) == 1);

  const ElementId<2> e_l0x0_l1y0(block_id, make_array(l0i0, l1i0), grid_index);
  const ElementId<2> e_l0x0_l1y1(block_id, make_array(l0i0, l1i1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l1y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l1y1) == 1);

  const ElementId<2> e_l1x0_l1y0(block_id, make_array(l1i0, l1i0), grid_index);
  const ElementId<2> e_l1x1_l1y0(block_id, make_array(l1i1, l1i0), grid_index);
  const ElementId<2> e_l1x0_l1y1(block_id, make_array(l1i0, l1i1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l1y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l1y0) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l1y1) == 2);

  const ElementId<2> e_l2x0_l0y0(block_id, make_array(l2i0, l0i0), grid_index);
  const ElementId<2> e_l2x1_l0y0(block_id, make_array(l2i1, l0i0), grid_index);
  const ElementId<2> e_l2x2_l0y0(block_id, make_array(l2i2, l0i0), grid_index);
  const ElementId<2> e_l2x3_l0y0(block_id, make_array(l2i3, l0i0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l0y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l0y0) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l0y0) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l0y0) == 3);

  const ElementId<2> e_l2x0_l1y0(block_id, make_array(l2i0, l1i0), grid_index);
  const ElementId<2> e_l2x1_l1y0(block_id, make_array(l2i1, l1i0), grid_index);
  const ElementId<2> e_l2x2_l1y0(block_id, make_array(l2i2, l1i0), grid_index);
  const ElementId<2> e_l2x3_l1y0(block_id, make_array(l2i3, l1i0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l1y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l1y0) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l1y0) == 4);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l1y0) == 6);

  const ElementId<2> e_l2x0_l1y1(block_id, make_array(l2i0, l1i1), grid_index);
  const ElementId<2> e_l2x1_l1y1(block_id, make_array(l2i1, l1i1), grid_index);
  const ElementId<2> e_l2x2_l1y1(block_id, make_array(l2i2, l1i1), grid_index);
  const ElementId<2> e_l2x3_l1y1(block_id, make_array(l2i3, l1i1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l1y1) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l1y1) == 3);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l1y1) == 5);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l1y1) == 7);

  const ElementId<2> e_l0x0_l2y0(block_id, make_array(l0i0, l2i0), grid_index);
  const ElementId<2> e_l0x0_l2y1(block_id, make_array(l0i0, l2i1), grid_index);
  const ElementId<2> e_l0x0_l2y2(block_id, make_array(l0i0, l2i2), grid_index);
  const ElementId<2> e_l0x0_l2y3(block_id, make_array(l0i0, l2i3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l2y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l2y1) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l2y2) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l0x0_l2y3) == 3);

  const ElementId<2> e_l1x0_l2y0(block_id, make_array(l1i0, l2i0), grid_index);
  const ElementId<2> e_l1x0_l2y1(block_id, make_array(l1i0, l2i1), grid_index);
  const ElementId<2> e_l1x0_l2y2(block_id, make_array(l1i0, l2i2), grid_index);
  const ElementId<2> e_l1x0_l2y3(block_id, make_array(l1i0, l2i3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l2y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l2y1) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l2y2) == 4);
  CHECK(domain::z_curve_index_from_element_id(e_l1x0_l2y3) == 6);

  const ElementId<2> e_l1x1_l2y0(block_id, make_array(l1i1, l2i0), grid_index);
  const ElementId<2> e_l1x1_l2y1(block_id, make_array(l1i1, l2i1), grid_index);
  const ElementId<2> e_l1x1_l2y2(block_id, make_array(l1i1, l2i2), grid_index);
  const ElementId<2> e_l1x1_l2y3(block_id, make_array(l1i1, l2i3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l2y0) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l2y1) == 3);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l2y2) == 5);
  CHECK(domain::z_curve_index_from_element_id(e_l1x1_l2y3) == 7);

  const ElementId<2> e_l2x0_l2y0(block_id, make_array(l2i0, l2i0), grid_index);
  const ElementId<2> e_l2x1_l2y0(block_id, make_array(l2i1, l2i0), grid_index);
  const ElementId<2> e_l2x2_l2y0(block_id, make_array(l2i2, l2i0), grid_index);
  const ElementId<2> e_l2x3_l2y0(block_id, make_array(l2i3, l2i0), grid_index);
  const ElementId<2> e_l2x0_l2y1(block_id, make_array(l2i0, l2i1), grid_index);
  const ElementId<2> e_l2x1_l2y1(block_id, make_array(l2i1, l2i1), grid_index);
  const ElementId<2> e_l2x2_l2y1(block_id, make_array(l2i2, l2i1), grid_index);
  const ElementId<2> e_l2x3_l2y1(block_id, make_array(l2i3, l2i1), grid_index);
  const ElementId<2> e_l2x0_l2y2(block_id, make_array(l2i0, l2i2), grid_index);
  const ElementId<2> e_l2x1_l2y2(block_id, make_array(l2i1, l2i2), grid_index);
  const ElementId<2> e_l2x2_l2y2(block_id, make_array(l2i2, l2i2), grid_index);
  const ElementId<2> e_l2x3_l2y2(block_id, make_array(l2i3, l2i2), grid_index);
  const ElementId<2> e_l2x0_l2y3(block_id, make_array(l2i0, l2i3), grid_index);
  const ElementId<2> e_l2x1_l2y3(block_id, make_array(l2i1, l2i3), grid_index);
  const ElementId<2> e_l2x2_l2y3(block_id, make_array(l2i2, l2i3), grid_index);
  const ElementId<2> e_l2x3_l2y3(block_id, make_array(l2i3, l2i3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l2y0) == 0);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l2y0) == 1);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l2y0) == 4);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l2y0) == 5);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l2y1) == 2);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l2y1) == 3);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l2y1) == 6);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l2y1) == 7);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l2y2) == 8);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l2y2) == 9);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l2y2) == 12);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l2y2) == 13);
  CHECK(domain::z_curve_index_from_element_id(e_l2x0_l2y3) == 10);
  CHECK(domain::z_curve_index_from_element_id(e_l2x1_l2y3) == 11);
  CHECK(domain::z_curve_index_from_element_id(e_l2x2_l2y3) == 14);
  CHECK(domain::z_curve_index_from_element_id(e_l2x3_l2y3) == 15);
}

// Test the Z-curve index computed for ElementIds in 3D. Tests
// domain::z_curve_index_from_element_id.
void test_z_curve_index_from_element_id_3d() {
  // The Z-curve does not depend on the block ID or grid index, so the choice of
  // these two values for this test is arbitrary
  const size_t block_id = 0;
  const size_t grid_index = 0;

  // Create SegmentIds for 3D test case. Since the 1D and 2D test cases are
  // exhaustive and exhaustively testing 3D in the same way is a lot to hard
  // code, we instead test one interesting 3D case. Having the highest
  // refinement (z ref = 4) at least 2 levels higher than the next highest
  // refinement (x ref = 2) is an important test case where the highest bits
  // of the Z-curve index come from one dimension, not an interleaving of more
  // than one dimension's segment index bits. More concretely, the bits for the
  // Z-curve index will be z3 z2 z1 x1 z0 x0, where the z3 z2 bits are not
  // interleaved with another dimension's segment index.
  const size_t x_refinement = 2;
  const size_t y_refinement = 0;
  const size_t z_refinement = 4;

  SegmentId x0(x_refinement, 0);
  SegmentId x1(x_refinement, 1);
  SegmentId x2(x_refinement, 2);
  SegmentId x3(x_refinement, 3);

  SegmentId y0(y_refinement, 0);

  SegmentId z0(z_refinement, 0);
  SegmentId z1(z_refinement, 1);
  SegmentId z2(z_refinement, 2);
  SegmentId z3(z_refinement, 3);
  SegmentId z4(z_refinement, 4);
  SegmentId z5(z_refinement, 5);
  SegmentId z6(z_refinement, 6);
  SegmentId z7(z_refinement, 7);
  SegmentId z8(z_refinement, 8);
  SegmentId z9(z_refinement, 9);
  SegmentId z10(z_refinement, 10);
  SegmentId z11(z_refinement, 11);
  SegmentId z12(z_refinement, 12);
  SegmentId z13(z_refinement, 13);
  SegmentId z14(z_refinement, 14);
  SegmentId z15(z_refinement, 15);

  // Create ElementIds with 3D SegmentIds from refinement 0 to 2 and check that
  // ElementIds are mapped to their correct Z-curve index
  // e_000 denotes segment indices x=0, y=0, z=0
  const ElementId<3> e_000(block_id, make_array(x0, y0, z0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_000) == 0);
  const ElementId<3> e_100(block_id, make_array(x1, y0, z0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_100) == 1);
  const ElementId<3> e_200(block_id, make_array(x2, y0, z0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_200) == 4);
  const ElementId<3> e_300(block_id, make_array(x3, y0, z0), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_300) == 5);

  const ElementId<3> e_001(block_id, make_array(x0, y0, z1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_001) == 2);
  const ElementId<3> e_101(block_id, make_array(x1, y0, z1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_101) == 3);
  const ElementId<3> e_201(block_id, make_array(x2, y0, z1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_201) == 6);
  const ElementId<3> e_301(block_id, make_array(x3, y0, z1), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_301) == 7);

  const ElementId<3> e_002(block_id, make_array(x0, y0, z2), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_002) == 8);
  const ElementId<3> e_102(block_id, make_array(x1, y0, z2), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_102) == 9);
  const ElementId<3> e_202(block_id, make_array(x2, y0, z2), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_202) == 12);
  const ElementId<3> e_302(block_id, make_array(x3, y0, z2), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_302) == 13);

  const ElementId<3> e_003(block_id, make_array(x0, y0, z3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_003) == 10);
  const ElementId<3> e_103(block_id, make_array(x1, y0, z3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_103) == 11);
  const ElementId<3> e_203(block_id, make_array(x2, y0, z3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_203) == 14);
  const ElementId<3> e_303(block_id, make_array(x3, y0, z3), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_303) == 15);

  const ElementId<3> e_004(block_id, make_array(x0, y0, z4), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_004) == 16);
  const ElementId<3> e_104(block_id, make_array(x1, y0, z4), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_104) == 17);
  const ElementId<3> e_204(block_id, make_array(x2, y0, z4), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_204) == 20);
  const ElementId<3> e_304(block_id, make_array(x3, y0, z4), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_304) == 21);

  const ElementId<3> e_005(block_id, make_array(x0, y0, z5), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_005) == 18);
  const ElementId<3> e_105(block_id, make_array(x1, y0, z5), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_105) == 19);
  const ElementId<3> e_205(block_id, make_array(x2, y0, z5), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_205) == 22);
  const ElementId<3> e_305(block_id, make_array(x3, y0, z5), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_305) == 23);

  const ElementId<3> e_006(block_id, make_array(x0, y0, z6), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_006) == 24);
  const ElementId<3> e_106(block_id, make_array(x1, y0, z6), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_106) == 25);
  const ElementId<3> e_206(block_id, make_array(x2, y0, z6), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_206) == 28);
  const ElementId<3> e_306(block_id, make_array(x3, y0, z6), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_306) == 29);

  const ElementId<3> e_007(block_id, make_array(x0, y0, z7), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_007) == 26);
  const ElementId<3> e_107(block_id, make_array(x1, y0, z7), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_107) == 27);
  const ElementId<3> e_207(block_id, make_array(x2, y0, z7), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_207) == 30);
  const ElementId<3> e_307(block_id, make_array(x3, y0, z7), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_307) == 31);

  const ElementId<3> e_008(block_id, make_array(x0, y0, z8), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_008) == 32);
  const ElementId<3> e_108(block_id, make_array(x1, y0, z8), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_108) == 33);
  const ElementId<3> e_208(block_id, make_array(x2, y0, z8), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_208) == 36);
  const ElementId<3> e_308(block_id, make_array(x3, y0, z8), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_308) == 37);

  const ElementId<3> e_009(block_id, make_array(x0, y0, z9), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_009) == 34);
  const ElementId<3> e_109(block_id, make_array(x1, y0, z9), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_109) == 35);
  const ElementId<3> e_209(block_id, make_array(x2, y0, z9), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_209) == 38);
  const ElementId<3> e_309(block_id, make_array(x3, y0, z9), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_309) == 39);

  const ElementId<3> e_0010(block_id, make_array(x0, y0, z10), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0010) == 40);
  const ElementId<3> e_1010(block_id, make_array(x1, y0, z10), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1010) == 41);
  const ElementId<3> e_2010(block_id, make_array(x2, y0, z10), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2010) == 44);
  const ElementId<3> e_3010(block_id, make_array(x3, y0, z10), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3010) == 45);

  const ElementId<3> e_0011(block_id, make_array(x0, y0, z11), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0011) == 42);
  const ElementId<3> e_1011(block_id, make_array(x1, y0, z11), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1011) == 43);
  const ElementId<3> e_2011(block_id, make_array(x2, y0, z11), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2011) == 46);
  const ElementId<3> e_3011(block_id, make_array(x3, y0, z11), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3011) == 47);

  const ElementId<3> e_0012(block_id, make_array(x0, y0, z12), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0012) == 48);
  const ElementId<3> e_1012(block_id, make_array(x1, y0, z12), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1012) == 49);
  const ElementId<3> e_2012(block_id, make_array(x2, y0, z12), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2012) == 52);
  const ElementId<3> e_3012(block_id, make_array(x3, y0, z12), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3012) == 53);

  const ElementId<3> e_0013(block_id, make_array(x0, y0, z13), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0013) == 50);
  const ElementId<3> e_1013(block_id, make_array(x1, y0, z13), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1013) == 51);
  const ElementId<3> e_2013(block_id, make_array(x2, y0, z13), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2013) == 54);
  const ElementId<3> e_3013(block_id, make_array(x3, y0, z13), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3013) == 55);

  const ElementId<3> e_0014(block_id, make_array(x0, y0, z14), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0014) == 56);
  const ElementId<3> e_1014(block_id, make_array(x1, y0, z14), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1014) == 57);
  const ElementId<3> e_2014(block_id, make_array(x2, y0, z14), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2014) == 60);
  const ElementId<3> e_3014(block_id, make_array(x3, y0, z14), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3014) == 61);

  const ElementId<3> e_0015(block_id, make_array(x0, y0, z15), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_0015) == 58);
  const ElementId<3> e_1015(block_id, make_array(x1, y0, z15), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_1015) == 59);
  const ElementId<3> e_2015(block_id, make_array(x2, y0, z15), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_2015) == 62);
  const ElementId<3> e_3015(block_id, make_array(x3, y0, z15), grid_index);
  CHECK(domain::z_curve_index_from_element_id(e_3015) == 63);
}

// Test the segment indices computed from a Z-curve index and creating initial
// ElementIds in Z-curve index order. Tests
// domain::segment_indices_from_z_curve_index and
// domain::initial_element_ids_in_z_curve_order.
template <size_t Dim>
void test_z_curve_index_mapping(
    const std::array<size_t, Dim>& block_refinement_levels) {
  // The Z-curve does not depend on the block ID or grid index, so the choice of
  // these two values for this test is arbitrary
  const size_t block_id = 1;
  const size_t grid_index = 2;

  std::vector<ElementId<Dim>> element_ids_in_default_order =
      initial_element_ids(block_id, block_refinement_levels, grid_index);

  const size_t num_elements = element_ids_in_default_order.size();
  // Whether or not we have encountered certain result Z-curve indices
  // when transforming ElementIds of a block to their Z-curve indices
  std::vector<bool> z_curve_index_hit(num_elements);
  std::fill(z_curve_index_hit.begin(), z_curve_index_hit.end(), false);

  // Check that computing the Z-curve index from an ElementId and then
  // Segment indices from the computed Z-curve index matches the Segment indices
  // of the original ElementId (i.e. transforming to and back from Z-curve index
  // recovers the original Segment indices of an ElementId). This checks that
  // domain::z_curve_index_from_element_i and
  // domain::segment_indices_from_z_curve_index are consistent with each other
  // in that they are inverses.
  for (const auto& element_id : element_ids_in_default_order) {
    const size_t result_z_curve_index =
        domain::z_curve_index_from_element_id(element_id);
    const std::array<size_t, Dim> result_segment_indices =
        domain::segment_indices_from_z_curve_index(result_z_curve_index,
                                                   block_refinement_levels);

    std::array<size_t, Dim> expected_segment_indices;
    for (size_t i = 0; i < Dim; ++i) {
      gsl::at(expected_segment_indices, i) = element_id.segment_id(i).index();
    }

    CHECK(result_segment_indices == expected_segment_indices);

    z_curve_index_hit[result_z_curve_index] = true;
  }

  // Check that there is a 1:1 mapping of ElementId to Z-curve index
  for (const size_t index_hit : z_curve_index_hit) {
    CHECK(index_hit);
  }

  const std::vector<ElementId<Dim>> element_ids_in_z_curve_order =
      domain::initial_element_ids_in_z_curve_order(
          block_id, block_refinement_levels, grid_index);

  // Check that there is a 1:1 mapping of ElementIds in the original input list
  // to the output list of ElementIds in Z-curve order (i.e. we preserved the
  // input list of ElementIds)
  CHECK(element_ids_in_z_curve_order.size() == num_elements);
  for (size_t i = 0; i < num_elements; i++) {
    const ElementId<Dim> element_id_from_z_curve_order =
        element_ids_in_z_curve_order[i];
    bool element_id_found = false;
    for (size_t j = 0; j < num_elements; j++) {
      const ElementId<Dim> element_id_from_default_order =
          element_ids_in_default_order[j];
      if (element_id_from_z_curve_order == element_id_from_default_order) {
        element_id_found = true;
        break;
      }
    }
    CHECK(element_id_found);
  }

  // Check that the ElementIds are ordered by Z-curve index
  for (size_t i = 0; i < element_ids_in_z_curve_order.size(); i++) {
    const auto& element_id = element_ids_in_z_curve_order[i];
    CHECK(domain::z_curve_index_from_element_id(element_id) == i);
  }
}

// Test the segment indices computed from a Z-curve index and creating initial
// ElementIds in Z-curve index order for each refinement level permutation for
// dimensions 1, 2, and 3
void test_z_curve_index_mapping_for_refinement_levels(
    const size_t min_refinement_level, const size_t max_refinement_level) {
  std::array<size_t, 1> block_refinement_levels_1d{};
  std::array<size_t, 2> block_refinement_levels_2d{};
  std::array<size_t, 3> block_refinement_levels_3d{};

  for (size_t i = min_refinement_level; i < max_refinement_level + 1; i++) {
    gsl::at(block_refinement_levels_1d, 0) = i;
    gsl::at(block_refinement_levels_2d, 0) = i;
    gsl::at(block_refinement_levels_3d, 0) = i;

    test_z_curve_index_mapping(block_refinement_levels_1d);

    for (size_t j = min_refinement_level; j < max_refinement_level + 1; j++) {
      gsl::at(block_refinement_levels_2d, 1) = j;
      gsl::at(block_refinement_levels_3d, 1) = j;

      test_z_curve_index_mapping(block_refinement_levels_2d);

      for (size_t k = min_refinement_level; k < max_refinement_level + 1; k++) {
        gsl::at(block_refinement_levels_3d, 2) = k;

        test_z_curve_index_mapping(block_refinement_levels_3d);
      }
    }
  }
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.ZCurve", "[Domain][Unit]") {
  // Test transforming ElementId to Z-curve index
  test_z_curve_index_from_element_id_1d();
  test_z_curve_index_from_element_id_2d();
  test_z_curve_index_from_element_id_3d();

  const size_t min_refinement_level = 0;
  const size_t max_refinement_level = 3;
  // Test transforming Z-curve index to ElementId and creating initial
  // ElementIds in Z-curve index order
  test_z_curve_index_mapping_for_refinement_levels(min_refinement_level,
                                                   max_refinement_level);
}
