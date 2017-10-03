#include <math.h>

#ifndef MM_VECKY
#define MM_VECKY

#ifndef __CUDACC__
  #define __host__
  #define __device__
  #define __global__
#endif

namespace pyscan {

#define OPERATION_STRUCT(OP_NAME, OPERATION, RET_T)		 \
template <typename E1, typename E2, typename T, unsigned int N>  \
 struct OP_NAME {							\
  template<unsigned int I>						\
  __host__ __device__ inline static RET_T				\
  apply(E1 const& v1, E2 const& v2){					\
    return v1.template get<I>() OPERATION v2.template get<I>();	\
  }									\
};									\
template <typename E1, typename T, unsigned int N>			\
struct OP_NAME<E1, T, T, N> {						\
  template<unsigned int I>						\
  __host__ __device__ inline static RET_T				\
  apply(E1 const& v1, T const& v2){					\
    return v1.template get<I>() OPERATION v2;				\
  }									\
};									\
 template <typename E2, typename T, unsigned int N>			\
struct OP_NAME<T, E2, T, N> {						\
   template<unsigned int I>						\
   __host__ __device__ inline static RET_T				\
   apply(T const& v1, E2 const& v2){					\
     return v1 OPERATION v2.template get<I>();	\
  }									\
};

//BINARY_FUNCTION_STRUCT(VecExponent, pow, T)
OPERATION_STRUCT(OprAdd, +, T)
OPERATION_STRUCT(OprMult, *, T)
OPERATION_STRUCT(OprDiv, /, T)
OPERATION_STRUCT(OprDiff, -, T)
OPERATION_STRUCT(OprLess, <, bool)
OPERATION_STRUCT(OprGreater, >, bool)
OPERATION_STRUCT(OprEqual, ==, bool)
OPERATION_STRUCT(OprNEq, ==, bool)
OPERATION_STRUCT(OprGEq, >=, bool)
OPERATION_STRUCT(OprLEq, <=, bool)


//A dummy struct
struct EmptyExpr {};


//This is for binary expressions
template <typename E1, template <typename, typename, typename, unsigned int> class OP, typename E2, typename T, unsigned int N>
struct ExprNode {
  E1 const& _e1;
  E2 const& _e2;
  __host__ __device__ inline ExprNode(E1 const& e1, E2 const& e2)
    : _e1(e1), _e2(e2)
  {}

  template<unsigned int I>
  __host__ __device__ inline T get() const {
    return OP<E1, E2, T, N>::template apply<I>(_e1, _e2);
  }
};

//This is for singleton expressions like negation or logical operations


template <typename E1, template <typename, typename, typename, unsigned int> class OP, typename T, unsigned int N>
struct ExprNode<E1, OP, EmptyExpr, T, N> {
  E1 const& _e1;
  __host__ __device__ inline ExprNode(E1 const& e1)
    : _e1(e1)
  {}

  template<unsigned int I>
  __host__ __device__ inline T get() const {
    return OP<E1, EmptyExpr, T, N>::template apply<I>(_e1);
  }
};

template<typename E1,
	 template <typename, typename, typename, unsigned int> class OP,
	 typename E2,
	 typename T,
	 unsigned int I,
	 unsigned int N>
struct COPY {
  __host__ __device__ inline static void
  apply(T* _data, ExprNode<E1, OP, E2, T, N> const& expr){
    _data[I] = expr.template get<I>();
    COPY<E1, OP, E2, T, I + 1, N>::apply(_data, expr);
  }
};

template<typename E1,
	 template <typename, typename, typename, unsigned int> class OP,
	 typename E2,
	 typename T,
	 unsigned int N>
struct COPY<E1, OP, E2, T, N, N> {
  __host__ __device__ inline static void
  apply(T* _data, ExprNode<E1, OP, E2, T, N> const& expr){
    (void)_data;
    (void)expr;
  }
};


template <typename T, unsigned int N>
class VecN {
  T _data[N];

  __host__ __device__
  template<typename Coord>
  void fill(int i, Coord f_el)  {
    _data[i] = f_el;
  }

  __host__ __device__
  template<typename Coord, typename... Coords>
  void fill(int i, Coord f_el, Coords... args) {
    _data[i] = f_el;
    fill(i - 1, args...);
  }

public:
  __host__ __device__
  template<typename... Coords>
  VecN(Coords... args) {
    fill(N - 1, args...);
  }

  template<typename E1, template <typename, typename, typename, unsigned int> class OP, typename E2>
  __host__ __device__ inline
  VecN(ExprNode<E1, OP, E2, T, N> const& expr){
    COPY<E1, OP, E2, T, 0, N>::apply(&_data[0], expr);
  }

  __host__ __device__ inline
  VecN(){}

  __host__ __device__ inline
  VecN(T const& val){
    for (unsigned int i = 0; i < N; i++){
      _data[i] = val;
    }
  }

  __host__ __device__ inline T& operator[](unsigned int i){
    return _data[i];
  }

  __host__ __device__ inline T operator[](unsigned int i) const {
    return _data[i];
  }

  template<unsigned int i>
  __host__ __device__ inline T get() const {
    return _data[i];
  }

  template<typename Expr>
  __host__ __device__ inline VecN<T, N>&
  operator += (Expr const& val){
    *this = *this + val;
    return *this;
  }

  template<typename Expr>
  __host__ __device__ inline VecN<T, N>&
  operator -= (Expr const& val){
    *this = *this - val;
    return *this;
  }

  template<typename Expr>
  __host__ __device__ inline VecN<T, N>&
  operator *= (Expr const& val){
    *this = *this * val;
    return *this;
  }

  template<typename Expr>
  __host__ __device__ inline VecN<T, N>&
  operator /= (Expr const& val){
    *this = *this / val;
    return *this;
  }
};

/////////////////////////////////////////////////////////////////////////////////////
/*
  This is for Unary vector operations that take a vector and apply a function to every element
 */
/////////////////////////////////////////////////////////////////////////////////////
#define UNARY_OPERATION_STRUCT(OP_NAME, OPERATION, RET_T)		\
template <typename E1, typename E2, typename T, unsigned int N>	        \
 struct OP_NAME {							\
  template<unsigned int I>						\
  __host__ __device__ inline static RET_T				\
  apply(E1 const& v1){							\
    return OPERATION v1.template get<I>();				\
  }									\
};


#define UNARY_FUNCTION_STRUCT(OP_NAME, FUNCTION, RET_T)		\
  template <typename E1, typename E2, typename T, unsigned int N>			\
 struct OP_NAME {							\
  template<unsigned int I>						\
  __host__ __device__ inline static RET_T				\
  apply(E1 const& v1){					\
    return (RET_T)FUNCTION(v1.template get<I>());				\
  }									\
};


#define UNARY_VECTOR_OPERATION(OP_NAME, SYMBOL, RET_T, F_NAME)		\
template<typename E1,							\
         typename E2,							\
	 template <typename, typename, typename, unsigned int> class PrevOp1, \
	 typename T,    \
	 unsigned int N>						\
__host__ __device__ inline ExprNode<ExprNode<E1, PrevOp1, E2, T, N>, OP_NAME, EmptyExpr, RET_T, N> const \
F_NAME SYMBOL (ExprNode<E1, PrevOp1, E2, T, N> const& expNode1){ \
  return ExprNode<ExprNode<E1, PrevOp1, E2, T, N>, OP_NAME, EmptyExpr, RET_T, N>(expNode1); \
}   \
template<typename T, unsigned int N>  \
__host__ __device__ inline ExprNode<VecN<T, N> , OP_NAME, EmptyExpr, RET_T, N> const  \
F_NAME SYMBOL (VecN<T, N> const& vec1){   \
  return ExprNode<VecN<T, N>, OP_NAME, EmptyExpr, RET_T, N>(vec1);  \
}


UNARY_OPERATION_STRUCT(UNegate, -, T)
UNARY_OPERATION_STRUCT(UBitNot, ~, T)
UNARY_OPERATION_STRUCT(ULogNot, !, bool)
UNARY_VECTOR_OPERATION(UNegate, -, T, operator)
UNARY_VECTOR_OPERATION(UBitNot, ~, T, operator)
UNARY_VECTOR_OPERATION(ULogNot, !, bool, operator)

#define UNARY_VECTOR_FUNCTION(F_NAME, RET_T)		\
  UNARY_FUNCTION_STRUCT(U##F_NAME, std::F_NAME, RET_T) \
  UNARY_VECTOR_OPERATION(U##F_NAME, std::F_NAME, RET_T, )
/*
UNARY_VECTOR_FUNCTION(cos, T)
UNARY_VECTOR_FUNCTION(sin, T)
UNARY_VECTOR_FUNCTION(tan, T)
UNARY_VECTOR_FUNCTION(acos, T)
UNARY_VECTOR_FUNCTION(asin, T)
UNARY_VECTOR_FUNCTION(atan, T)
UNARY_VECTOR_FUNCTION(cosh, T)
UNARY_VECTOR_FUNCTION(sinh, T)
UNARY_VECTOR_FUNCTION(tanh, T)
UNARY_VECTOR_FUNCTION(log, T)
UNARY_VECTOR_FUNCTION(log10, T)
UNARY_VECTOR_FUNCTION(sqrt, T)
UNARY_VECTOR_FUNCTION(ceil, T)
UNARY_VECTOR_FUNCTION(floor, T)
UNARY_VECTOR_FUNCTION(fmod, T)

UNARY_VECTOR_FUNCTION(abs, T)
*/



/////////////////////////////////////////////////////////////////////////////////////
/*
  This is for vector/vector operations that produce vectors.
 */
/////////////////////////////////////////////////////////////////////////////////////

#define BINARY_VECTOR_OPERATION(OP_NAME, SYMBOL, RET_T, F_NAME)		\
template<typename E1,							\
         typename E2,							\
	 typename E3,                                                   \
	 typename E4,  \
	 template <typename, typename, typename, unsigned int> class PrevOp1, \
	 template <typename, typename, typename, unsigned int> class PrevOp2,  \
	 typename T,    \
	 unsigned int N>						\
__host__ __device__ inline ExprNode<ExprNode<E1, PrevOp1, E2, T, N>, OP_NAME, ExprNode<E3, PrevOp2, E4, T, N>, RET_T, N> const \
F_NAME SYMBOL (ExprNode<E1, PrevOp1, E2, T, N> const& expNode1, ExprNode<E3, PrevOp2, E4, T, N> const& expNode2){ \
  return ExprNode<ExprNode<E1, PrevOp1, E2, T, N>, OP_NAME, ExprNode<E3, PrevOp2, E4, T, N>, RET_T, N>(expNode1, expNode2); \
}   \
    \
template<typename E1,   \
	 typename E2,   \
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,  \
	 unsigned int N>						\
__host__ __device__ inline ExprNode<ExprNode<E1, PrevOp, E2, T, N>, OP_NAME, VecN<T, N>, RET_T, N> const   \
F_NAME SYMBOL (ExprNode<E1, PrevOp, E2, T, N> const& expNode, VecN<T, N> const& vec){  \
  return ExprNode<ExprNode<E1, PrevOp, E2, T, N>, OP_NAME, VecN<T, N>, RET_T, N>(expNode, vec);  \
}  \
  \
template<typename E1,  \
	 typename E2, \
	 template <typename, typename, typename, unsigned int> class PrevOp,  \
	 typename T,   \
	 unsigned int N>						\
__host__ __device__ inline ExprNode<VecN<T, N>, OP_NAME, ExprNode<E1, PrevOp, E2, T, N>, RET_T, N> const  \
F_NAME SYMBOL (VecN<T, N> const& vec, ExprNode<E1, PrevOp, E2, T, N> const& expNode){  \
  return ExprNode<VecN<T, N>, OP_NAME, ExprNode<E1, PrevOp, E2, T, N>, RET_T, N>(vec, expNode);  \
}  \
  \
template<typename T, unsigned int N>  \
__host__ __device__ inline ExprNode<VecN<T, N> , OP_NAME, VecN<T, N>, RET_T, N> const  \
F_NAME SYMBOL (VecN<T, N> const& vec1, VecN<T, N> const& vec2){   \
  return ExprNode<VecN<T, N>, OP_NAME, VecN<T, N>, RET_T, N>(vec1, vec2);  \
}


BINARY_VECTOR_OPERATION(OprAdd, +, T, operator)
BINARY_VECTOR_OPERATION(OprMult, *, T, operator)
BINARY_VECTOR_OPERATION(OprDiv, /, T, operator)
BINARY_VECTOR_OPERATION(OprDiff, -, T, operator)
BINARY_VECTOR_OPERATION(OprLess, <, bool, operator)
BINARY_VECTOR_OPERATION(OprGreater, >, bool, operator)
BINARY_VECTOR_OPERATION(OprEqual, ==, bool, operator)
BINARY_VECTOR_OPERATION(OprNEq, !=, bool, operator)
BINARY_VECTOR_OPERATION(OprGEq, ==, bool, operator)
BINARY_VECTOR_OPERATION(OprLEq, !=, bool, operator)

template <typename E1, typename E2, typename T, unsigned int N>
struct CROSS {
  template <unsigned int I>
  __host__ __device__ inline static T
  apply(E1 const& v1, E2 const& v2){
    return (v1.template get<(I  + 1) % 3>() * v2.template get<(I + 2) % 3>() -
	    v1.template get<(I  + 2) % 3>() * v2.template get<(I + 1) % 3>());
  }
};

BINARY_VECTOR_OPERATION(CROSS, , T, cross)


/////////////////////////////////////////////////////////////////////////////////////
/*
  This is for vector/scalar operations that produce vectors.
 */
/////////////////////////////////////////////////////////////////////////////////////

#define SCALAR_VECTOR_OPERATION(OP_NAME, SYMBOL, RET_T, F_NAME)		\
template<typename E1,							\
	 typename E2,							\
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,							\
	 unsigned int N>						\
__host__ __device__ inline  ExprNode<ExprNode<E1, PrevOp, E2, T, N>, OP_NAME, T, RET_T, N> const \
F_NAME SYMBOL (ExprNode<E1, PrevOp, E2, T, N> const& expNode, T const& val){	\
  return ExprNode<ExprNode<E1, PrevOp, E2, T, N>, OP_NAME, T, T, N>(expNode, val); \
}    \
template<typename E1,   \
	 typename E2,     \
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,  \
	 unsigned int N>						\
__host__ __device__ inline ExprNode<T, OP_NAME, ExprNode<E1, PrevOp, E2, T, N>, RET_T, N> const \
F_NAME SYMBOL (T const& val, ExprNode<E1, PrevOp, E2, T, N> const& expNode){  \
  return ExprNode<T, OP_NAME, ExprNode<E1, PrevOp, E2, T, N>, T, N>(val, expNode); \
}									\
template<typename T,							\
	 unsigned int N>						\
__host__ __device__ inline  ExprNode<VecN<T, N>, OP_NAME, T, RET_T, N> const	\
F_NAME SYMBOL (VecN<T, N> const& vec, T const& val){			\
  return ExprNode<VecN<T, N>, OP_NAME, T, T, N>(vec, val);		\
}									\
template<typename T,							\
	 unsigned int N>						\
__host__ __device__ inline  ExprNode<T, OP_NAME, VecN<T, N>, RET_T, N> const	\
F_NAME SYMBOL (T const& val, VecN<T, N> const& vec){			\
  return ExprNode<T, OP_NAME, VecN<T, N>, T, N>(val, vec); \
}

SCALAR_VECTOR_OPERATION(OprMult, *, T, operator)
SCALAR_VECTOR_OPERATION(OprDiv, /, T, operator)
SCALAR_VECTOR_OPERATION(OprLess, <, bool, operator)
SCALAR_VECTOR_OPERATION(OprGreater, >, bool, operator)
SCALAR_VECTOR_OPERATION(OprEqual, ==, bool, operator)
SCALAR_VECTOR_OPERATION(OprNEq, !=, bool, operator)
SCALAR_VECTOR_OPERATION(OprGEq, ==, bool, operator)
SCALAR_VECTOR_OPERATION(OprLEq, !=, bool, operator)


/////////////////////////////////////////////////////////////////////////////////////
/*
  This is for operations like dot product that take two vectors and produce a scalar
  value of some type.
 */
/////////////////////////////////////////////////////////////////////////////////////


#define BINARY_VECTOR_TO_SCALAR(OP_NAME, SYMBOL, RET_T, F_NAME)		\
template<typename E1,							\
	 typename E2,							\
	 typename E3,							\
	 typename E4,							\
	 template <typename, typename, typename, unsigned int> class PrevOp1, \
	 template <typename, typename, typename, unsigned int> class PrevOp2, \
	 typename T,							\
	 unsigned int N>						\
__host__ __device__ inline RET_T const					\
F_NAME SYMBOL (ExprNode<E1, PrevOp1, E2, T, N> const& expNode1, ExprNode<E3, PrevOp2, E4, T, N> const& expNode2){   \
  return OP_NAME<ExprNode<E1, PrevOp1, E2, T, N>, ExprNode<E3, PrevOp2, E4, T, N>, T, 0, N>::apply(expNode1, expNode2);   \
}   \
    \
template<typename E1,							\
	 typename E2,							\
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,							\
	 unsigned int N>						\
__host__ __device__ inline RET_T const					\
F_NAME SYMBOL (ExprNode<E1, PrevOp, E2, T, N> const& expNode, VecN<T, N> const& vec){  \
  return OP_NAME<ExprNode<E1, PrevOp, E2, T, N>, VecN<T, N>, T, 0, N>::apply(expNode, vec);  \
}									\
									\
template<typename E1,							\
	 typename E2,							\
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,\
	 unsigned int N>  \
__host__ __device__ inline RET_T const  \
F_NAME SYMBOL (VecN<T, N> const& vec, ExprNode<E1, PrevOp, E2, T, N> const& expNode){ \
  return OP_NAME<VecN<T, N>, ExprNode<E1, PrevOp, E2, T, N>, T, 0, N>::apply(vec, expNode); \
}									\
 									\
template<typename T,							\
	 unsigned int N>						\
__host__ __device__ inline T const					\
F_NAME SYMBOL (VecN<T, N> const& vec1, VecN<T, N> const& vec2){	        \
  return OP_NAME<VecN<T, N>, VecN<T, N>, T, 0, N>::apply(vec1, vec2);	\
}

template<typename E1, typename E2, typename T, unsigned int I, unsigned int N>
struct DOT {
  __host__ __device__ inline static T
  apply(E1 const& e1, E2 const& e2) {
    return e1.template get<I>() * e2.template get<I>() + DOT<E1, E2, T, I + 1, N>::apply(e1, e2);
  }
};

template<typename E1, typename E2, typename T, unsigned int N>
struct DOT<E1, E2, T, N, N> {
  __host__ __device__ inline static T
  apply(E1 const& e1, E2 const& e2) {
    (void)e1;
    (void)e2;
    return 0.0;
  }
};

BINARY_VECTOR_TO_SCALAR(DOT, , T, dot)


#define UNARY_VECTOR_TO_SCALAR(OP_NAME, SYMBOL, RET_T, F_NAME)   \
template<typename E1,							\
	 typename E2,							\
	 template <typename, typename, typename, unsigned int> class PrevOp, \
	 typename T,							\
	 unsigned int N>						\
__host__ __device__ inline RET_T const					\
F_NAME SYMBOL (ExprNode<E1, PrevOp, E2, T, N> const& expNode){	        \
  return OP_NAME<ExprNode<E1, PrevOp, E2, T, N>, T, 0, N>::apply(expNode);   \
}   \
template<typename T,							\
	 unsigned int N>						\
__host__ __device__ inline T const					\
F_NAME SYMBOL (VecN<T, N> const& vec){					\
  return OP_NAME<VecN<T, N>, T, 0, N>::apply(vec);			\
}


template<typename E1, typename T, unsigned int I, unsigned int N>
struct MAG {
  __host__ __device__ inline static T
  apply(E1 const& e1) {
    return sqrt(DOT<E1, E1, T, I, N>::apply(e1, e1));
  }
};

template<typename E1, typename T, unsigned int I, unsigned int N>
struct SUM {
  __host__ __device__ inline static T
  apply(E1 const& e1) {
    return e1.template get<I>() + SUM<E1, T, I + 1, N>::apply(e1);
  }
};

template<typename E1, typename T, unsigned int N>
struct SUM<E1, T, N, N> {
  __host__ __device__ inline static T
  apply(E1 const& e1) {
    (void)e1;
    return 0.0;
  }
};



UNARY_VECTOR_TO_SCALAR(MAG, , T, mag)
UNARY_VECTOR_TO_SCALAR(SUM, , T, sum)


using VecD = VecN<double, 2>;

}
#endif
