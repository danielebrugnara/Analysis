#include <ReactionReconstruction.h>
//Force compilation of following templates, necessary to have implementation of
//class in file different from header
template class ReactionReconstruction<float>;
template class ReactionReconstruction2body<float>;
template class ReactionReconstruction3body<float>;

template class ReactionReconstruction<double>;
template class ReactionReconstruction2body<double>;
template class ReactionReconstruction3body<double>;

template class ReactionReconstruction<long double>;
template class ReactionReconstruction2body<long double>;
template class ReactionReconstruction3body<long double>;
