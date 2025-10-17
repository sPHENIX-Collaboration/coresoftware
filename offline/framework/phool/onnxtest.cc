#include <onnxruntime_cxx_api.h>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[])
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " model.onnx" << std::endl;
    return 1;
  }

  const std::string model_path = argv[1];

  try
  {
    Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "meta_reader");
    Ort::SessionOptions session_options;
    session_options.SetIntraOpNumThreads(1);

    // Create session to load the model
    Ort::Session session(env, model_path.c_str(), session_options);

    std::cout << "✅ Model loaded successfully: " << model_path << "\n";
    std::cout << "----------------------------------------\n";

    // Print model input information
    Ort::AllocatorWithDefaultOptions allocator;

    size_t num_input_nodes = session.GetInputCount();
    std::vector<std::string> inputnames = session.GetInputNames();
    std::cout << "Inputs (" << num_input_nodes << "):\n";
    for (size_t i = 0; i < num_input_nodes; i++)
    {
      auto type_info = session.GetInputTypeInfo(i);
      auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

      ONNXTensorElementDataType type = tensor_info.GetElementType();
      auto input_dims = tensor_info.GetShape();

      std::cout << "  • " << inputnames[i] << " (type=" << type << ", shape=[";
      for (size_t j = 0; j < input_dims.size(); j++)
      {
        std::cout << input_dims[j];
        if (j + 1 < input_dims.size())
        {
          std::cout << ", ";
        }
      }
      std::cout << "])\n";
    }

    // Print model output information
    size_t num_output_nodes = session.GetOutputCount();
    std::cout << "\nOutputs (" << num_output_nodes << "):\n";
    std::vector<std::string> outputnames = session.GetOutputNames();
    for (size_t i = 0; i < num_output_nodes; i++)
    {
      auto type_info = session.GetOutputTypeInfo(i);
      auto tensor_info = type_info.GetTensorTypeAndShapeInfo();

      ONNXTensorElementDataType type = tensor_info.GetElementType();
      auto output_dims = tensor_info.GetShape();

      std::cout << "  • " << outputnames[i] << " (type=" << type << ", shape=[";
      for (size_t j = 0; j < output_dims.size(); j++)
      {
        std::cout << output_dims[j];
        if (j + 1 < output_dims.size())
        {
          std::cout << ", ";
        }
      }
      std::cout << "])\n";
    }

    std::cout << "\n----------------------------------------\n";
    std::cout << "ONNX Runtime Version: " << Ort::GetVersionString() << "\n";
  }
  catch (const Ort::Exception& e)
  {
    std::cerr << "❌ ONNX Runtime Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
