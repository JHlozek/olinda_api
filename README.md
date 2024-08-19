# olinda_model_runner

### Load Olinda ONNX model

```
import onnx
import onnxruntime as rt
onnx_model = onnx.load("/path/to/model.onnx")
```

### Setup ONNX runtime

```
onnx_rt = rt.InferenceSession(onnx_model.SerializeToString())
output_names = [n.name for n in onnx_model.graph.output]

```

### Featurize SMILES and run prediction

```
from morgan_featurizer import MorganFeaturizer
featurizer = MorganFeaturizer()
fingerprint = featurizer.featurize(["CCC"])
onnx_rt.run(output_names, {"input": fingerprint})
```
