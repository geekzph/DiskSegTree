#### 简介

通过线段树查询最大字段和； 将原始的二叉线段树修改为类似B树的多分支树，存储在存盘上，以减少I/O次数；

#### 编译方法

使用Visual Studio 2013 编译通过；

编译时需要修改 Project Property -> C/C++ -> Precompile Headers 为 Not Using Precompiled Headers

#### 参数设置

kBranchNum为线段树的分支数
