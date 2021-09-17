# **4DFoot**
This is a repository for the paper "4D-Foot: A fully automated pipeline of four-dimensional analysis of the foot bones using bi-plane x-ray video and CT" accepted by MICCAI 2021. This repository depends on https://github.com/YoshitoOtake/RegTools.

![alt text](https://github.com/YoshitoOtake/4DFoot/tree/main/figs/Fig_MethodOverview_v2.png?raw=true "Method Overview")

## Abstract
We aim to elucidate the mechanism of the foot by automated measurement of its multiple bone movement using 2D-3D registration of bi-plane x-ray video and a stationary 3D CT. Conventional analyses allowed tracking of only 3 large proximal tarsal bones due to the requirement of manual segmentation and manual initialization of 2D-3D registration. The learning-based 2D-3D registration, on the other hand, has been actively studied and demonstrating a large capture range, but the accuracy is inferior to conventional optimization-based methods. We propose a fully automated pipeline using a cost function that seamlessly incorporates the reprojection error at the landmarks in CT and x-ray detected by off-the-shelf CNNs into the conventional image similarity cost, combined with the automated bone segmentation. We experimentally demonstrated that the pipeline allowed a robust and accurate 2D-3D registration to track all 12 tarsal bones, including the metatarsals at the foot arch, which is especially important in the foot biomechanics but has been unmeasurable with previous methods. We evaluated the proposed fully automated pipeline in studies using a bone phantom and real x-ray images of human subjects. The real image study showed the registration error of 0.38 ± 1.95 mm in translation and 0.38 ± 1.20 degrees in rotation for the proximal tarsal bones.

## Citation

```
@inproceedings{otake2021,
  title={4D-Foot: A fully automated pipeline of four-dimensional analysis of the foot bones using bi-plane x-ray video and CT},
  author={Shuntaro Mizoe, Yoshito Otake, Takuma Miyamoto, Mazen Soufi, Satoko Nakao, Yasuhito Tanaka, and Yoshinobu Sato},
  booktitle={International Conference on Medical Image Computing and Computer-Assisted Intervention},
  year={2021},
  organization={Springer}
}
```
