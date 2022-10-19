#ifndef ONNXLIB_H
#define ONNXLIB_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wpedantic"
#include <onnxruntime_cxx_api.h>
#pragma GCC diagnostic pop

// This is a stub for some ONNX code refactoring

Ort::Session *onnxSession(std::string &modelfile);

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nsamp, int Nreturn);

#endif
