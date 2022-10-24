#include "onnxlib.h"

#include <iostream>

// --------------------------------------------------
Ort::Session *onnxSession(std::string &modelfile)
{
  Ort::Env env(OrtLoggingLevel::ORT_LOGGING_LEVEL_WARNING, "fit");
  Ort::SessionOptions sessionOptions;
  sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);

  return new Ort::Session(env, modelfile.c_str(), sessionOptions);
}

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nsamp, int Nreturn)
{
  Ort::MemoryInfo memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

  Ort::AllocatorWithDefaultOptions allocator;

  std::vector<Ort::Value> inputTensors, outputTensors;

  std::vector<int64_t> inputDimsN = {N, Nsamp};
  std::vector<int64_t> outputDimsN = {N, Nreturn};

  std::vector<float> outputTensorValuesN(N * Nreturn);

  inputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, input.data(), N * Nsamp, inputDimsN.data(), inputDimsN.size()));
  outputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, outputTensorValuesN.data(), N * Nreturn, outputDimsN.data(), outputDimsN.size()));

  std::vector<const char *> inputNames{session->GetInputName(0, allocator)};
  std::vector<const char *> outputNames{session->GetOutputName(0, allocator)};

  session->Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);

  return outputTensorValuesN;
}
