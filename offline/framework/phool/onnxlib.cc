#include "onnxlib.h"

#include <iostream>

Ort::Session *onnxSession(std::string &modelfile, int verbosity)
{
  Ort::Env env(OrtLoggingLevel::ORT_LOGGING_LEVEL_WARNING, "fit");
  Ort::SessionOptions sessionOptions;
  sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_EXTENDED);
  auto *session = new Ort::Session(env, modelfile.c_str(), sessionOptions);
  auto type_info = session->GetInputTypeInfo(0);
  auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
  auto input_dims = tensor_info.GetShape();
  onnxlib::n_input = input_dims[1];
  type_info = session->GetOutputTypeInfo(0);
  tensor_info = type_info.GetTensorTypeAndShapeInfo();
  auto output_dims = tensor_info.GetShape();
  onnxlib::n_output = output_dims[1];
  if (verbosity > 0)
  {
    std::cout << "onnxlib: using model " << modelfile << std::endl;
    std::cout << "Number of Inputs: " << onnxlib::n_input << std::endl;
    std::cout << "Number of Outputs: " << onnxlib::n_output << std::endl;
  }
  return session;
}

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nsamp, int Nreturn)
{
  Ort::MemoryInfo memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

  Ort::AllocatorWithDefaultOptions allocator;

  std::vector<Ort::Value> inputTensors;
  std::vector<Ort::Value> outputTensors;

  std::vector<int64_t> inputDimsN = {N, Nsamp};
  std::vector<int64_t> outputDimsN = {N, Nreturn};
  int inputlen = N * Nsamp;
  int outputlen = N * Nreturn;

  std::vector<float> outputTensorValuesN(outputlen);

  inputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, input.data(), inputlen, inputDimsN.data(), inputDimsN.size()));
  outputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, outputTensorValuesN.data(), outputlen, outputDimsN.data(), outputDimsN.size()));

#if ORT_API_VERSION == 12
  std::vector<const char *> inputNames{session->GetInputName(0, allocator)};
  std::vector<const char *> outputNames{session->GetOutputName(0, allocator)};
#elif ORT_API_VERSION == 22
  std::vector<const char *> inputNames;
  std::vector<const char *> outputNames;

  char *name{nullptr};
  for (const std::string &s : session->GetInputNames())
  {
    name = new char[s.size() + 1];
    sprintf(name, "%s", s.c_str()); //NOLINT(hicpp-vararg)
    inputNames.push_back(name);
  }
  for (const std::string &s : session->GetOutputNames())
  {
    name = new char[s.size() + 1];
    sprintf(name, "%s", s.c_str()); //NOLINT(hicpp-vararg)
    outputNames.push_back(name);
  }
#else
#define XSTR(x) STR(x)
#define STR(x) #x
#pragma message "ORT_API_VERSION " XSTR(ORT_API_VERSION) " not implemented"
#endif
  session->Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);

#if ORT_API_VERSION == 22
  for (const auto *iter : inputNames)
  {
    delete[] iter;
  }
  for (const auto *iter : outputNames)
  {
    delete[] iter;
  }
#endif
  return outputTensorValuesN;
}

std::vector<float> onnxInference(Ort::Session *session, std::vector<float> &input, int N, int Nx, int Ny, int Nz, int Nreturn)
{
  // Define the memory information for ONNX Runtime
  Ort::MemoryInfo memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

  Ort::AllocatorWithDefaultOptions allocator;

  std::vector<int64_t> inputDims = {N, Nx, Ny, Nz};
  std::vector<int64_t> outputDimsN = {N, Nreturn};
  int inputlen = N * Nx * Ny * Nz;
  int outputlen = N * Nreturn;

  std::vector<float> outputTensorValues(outputlen);

  std::vector<Ort::Value> inputTensors;
  std::vector<Ort::Value> outputTensors;

  inputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, input.data(), inputlen, inputDims.data(), inputDims.size()));

  outputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, outputTensorValues.data(), outputlen, outputDimsN.data(), outputDimsN.size()));

#if ORT_API_VERSION == 12
  std::vector<const char *> inputNames{session->GetInputName(0, allocator)};
  std::vector<const char *> outputNames{session->GetOutputName(0, allocator)};
#elif ORT_API_VERSION == 22
  std::vector<const char *> inputNames;
  std::vector<const char *> outputNames;
  char *name{nullptr};
  for (const std::string &s : session->GetInputNames())
  {
    name = new char[s.size() + 1];
    sprintf(name, "%s", s.c_str()); //NOLINT(hicpp-vararg)
    inputNames.push_back(name);
  }
  for (const std::string &s : session->GetOutputNames())
  {
    name = new char[s.size() + 1];
    sprintf(name, "%s", s.c_str()); //NOLINT(hicpp-vararg)
    outputNames.push_back(name);
  }
#else
#define XSTR(x) STR(x)
#define STR(x) #x
#pragma message "ORT_API_VERSION " XSTR(ORT_API_VERSION) " not implemented"
#endif
  session->Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);
#if ORT_API_VERSION == 22
  for (const auto *iter : inputNames)
  {
    delete[] iter;
  }
  for (const auto *iter : outputNames)
  {
    delete[] iter;
  }
#endif

  return outputTensorValues;
}
