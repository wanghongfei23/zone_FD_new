
### 章节结构及内容简介
#### 1. INTRODUCTION
- **内容**：介绍CGNS（计算流体动力学通用标记系统）的起源、目标及组成。CGNS由波音和NASA于1994年联合开发，旨在标准化CFD数据的输入输出，支持结构化/非结构化网格、流动解、边界条件等数据存储，采用ADF/HDF5作为数据库管理器，提供底层CGIO接口和中层API库。此外，还介绍了CGNS文件的树状节点结构及本用户指南的组织方式。

#### 2. GETTING STARTED
- **2.1 Structured Grid**：通过实例演示结构化网格的创建与操作，包括单区网格、带流动解的单区网格、带边界条件的单区网格，以及带1-to-1连接的多区网格，涉及FORTRAN和C代码示例，说明网格坐标、流动解数据的读写方法。
- **2.2 Unstructured Grid**：介绍非结构化网格的处理，包括单区网格、带流动解的单区网格及带边界条件的单区网格，强调非结构化网格中元素连接和边界条件定义的特点。

#### 3. ADDITIONAL INFORMATION
- **3.1 Convergence History**：说明如何存储CFD解的收敛历史数据，如升力系数随迭代次数的变化，数据存储在CGNSBase层级的ConvergenceHistory节点。
- **3.2 Descriptor Nodes**：介绍描述节点的用途，可用于插入注释或存储ASCII输入文件，支持包含换行符的文本字符串。
- **3.3 Dimensional Data**：阐述维度数据的表示方法，通过DataClass、DimensionalUnits和DimensionalExponents节点定义数据单位和量纲，如MKS单位下的密度和压力。
- **3.4 Nondimensional Data**：区分两种无量纲数据类型（NormalizedByDimensional和NormalizedByUnknownDimensional），说明参考状态节点（ReferenceState）的作用，如定义马赫数和雷诺数。
- **3.5 Flow Equation Sets**：描述流动方程集的记录方式，如使用FlowEquationSet节点记录求解的方程类型（如纳维-斯托克斯方程）、湍流模型和气体模型。
- **3.6 Time-Dependent Data**：介绍时间相关数据的存储，通过BaseIterativeData和ZoneIterativeData节点关联多个流动解，支持非移动网格的时间精确模拟。
- **3.7 Using Links**：说明链接的使用，可在同一文件或不同文件间关联节点，避免重复存储数据，如多个流动解共享同一网格坐标。

#### 4. TROUBLESHOOTING
- **4.1 Handling Errors**：提及API的错误检查机制，建议用户检查错误代码并优雅退出程序。
- **4.2 Known Problems**：指出链接使用中的限制，如修改模式下链接文件需有写权限。

#### 5. FREQUENTLY ASKED QUESTIONS
- **内容**：解答常见问题，包括自定义数据名称的兼容性、Family的作用、DiscreteData和IntegralData的用途、编程最佳实践、CGNS文件查看方法及合规性检查方式。

#### 附录A. EXAMPLE COMPUTER CODES
- **内容**：列出可从SourceForge获取的示例代码，涵盖结构化网格、非结构化网格、收敛历史等场景的FORTRAN和C代码，如write_grid_str.f和read_grid_str.f。

#### 附录B. OVERVIEW OF THE SIDS
- **B.1 The Big Picture**：概述SIDS（标准接口数据结构）的层次结构，说明CGNS文件的根节点、CGNSBase节点及Zone节点的组成，强调节点标签和名称的规范。
- **B.2 Implementation at the Lower Levels**：讨论底层数据存储规范，如使用标准化名称（如Density）和UserDefinedData节点处理特殊数据。
- **B.3 Boundary Conditions**：介绍边界条件的层次结构，支持简化实现和完整SIDS合规实现，如BCType和BCDataSet的组合使用。
- **B.4 Zone Connectivity**：说明区域连接的类型（1-to-1、 patched、overset），通过ZoneGridConnectivity节点及其子节点（如GridConnectivity1to1_t）定义连接关系。
- **B.5 Structured Zone Example**：以二维结构化网格为例，展示SIDS合规的文件布局，包括网格坐标、流动解和连接信息的存储方式。

#### 附录C. GUIDELINE FOR PLOT3D VARIABLES
- **内容**：提供PLOT3D格式数据在CGNS中的存储指南，区分维度数据、已知归一化无量纲数据和未知归一化无量纲数据的处理方式，确保数据可移植性，如定义参考状态下的密度和音速。
