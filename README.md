# Systimator-DataFlow
A Simple Dataflow for Systimator on Tiny-Yolo, Yolo, AlexNet, CaffeNet and VGG16. <br>
Pick any of the following architecture and run Systimator.
# Tiny-Yolo Network Architecture
<pre>
% Layer   Filters   Size   Stride   IFM    IFM_Depth   OFM   OFM_Depth
    1       16       3       1      416        3       416       16
    2       32       3       1      208       16       208       32
    3       64       3       1      104       32       104       64
    4      128       3       1       52       64        52      128
    5      256       3       1       26      128        26      256
    6      512       3       1       13      256        13      512
    7     1024       3       1       13      512        13     1024
    8     1024       3       1       13     1024        13     1024
    9      515       1       1       13     1024        13      515
</pre>
# Yolo Network Architecture
<pre>
    %Filters   Size        FM    IFM_Depth
     32       3           416        3
     64       3           208       32
     128      3           104       64
     64       1           104      128
    128       3           104       64
    256       3            52      128
    128       1            52      256
    256       3            52      128
    512       3            26      256
    256       1            26      512
    512       3            26      256
    256       1            26      512
    512       3            26      256
   1024       3            13      512
    512       1            13     1024
   1024       3            13      512
    512       1            13     1024
   1024       3            13      512
   1024       3            13     1024
   1024       3            13     1024
   </pre>
## AlexNet Network Architecture
<pre>
  Layer   Filters   Size   Stride   IFM    IFM_Depth   OFM   OFM_Depth
    1       96      11       4      227        3        55       96
    2      256       5       1       27       96        27      256
    3      384       3       1       13      256        13      384
    4      384       3       1       13      384        13      384
    5      256       3       1       13      384        13      256
</pre>
# CaffeNet Network Architecture
<pre>
  Filters   Size        FM    IFM_Depth
    96       11         227        3
    256       5          55       96
    384       3          27      256
    384       3          13      384
    256       3          13      256
</pre>
# VGG16 Network Architecture
<pre>
    %Filters   Size        FM    IFM_Depth
     64         3          224     1
     64         3          224     64
    128         3          112     64
    128         3          112    128
    256         3           56    128
    256         3           56    256
    256         3           56    256
    512         3           28    256
    512         3           28    512
    512         3           28    512
    512         3           14    512
    512         3           14    512
    512         3           14    512
    </pre>


