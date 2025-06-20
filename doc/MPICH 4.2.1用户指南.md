
### 章节结构及内容简介  
#### 1. INTRODUCTION  
- **内容**：本手册假设MPICH已安装，主要介绍如何编译、链接和运行MPI应用程序，以及使用MPICH附带的工具。手册为初步版本，部分内容待完善，但已足够帮助用户入门。  

#### 2. Getting Started with MPICH  
- **2.1 Default Runtime Environment**：MPICH默认运行环境为Hydra，支持进程管理与通信分离，也可使用其他进程管理器。  
- **2.2 Starting Parallel Jobs**：MPICH实现了mpiexec及其标准参数，并提供针对不同进程管理系统的扩展，详见第5章。  
- **2.3 Command-Line Arguments in Fortran**：MPICH不强制Fortran程序访问命令行参数，若需要，需自行确保Fortran环境中iargc和getarg函数可用。  

#### 3. Quick Start  
- **内容**：指导用户将MPICH的bin目录添加到路径中，通过运行mpichversion查看版本，并使用mpiexec运行示例程序（如cpi），验证MPICH安装正确性。  

#### 4. Compiling and Linking  
- **4.1 Special Issues for C++**：C++编译时可能因stdio.h与MPI C++接口冲突导致错误，可通过#undef相关宏或添加-DMPICH_IGNORE_CXX_SEEK编译选项解决。  
- **4.2 Special Issues for Fortran**：Fortran 77使用mpif.h定义MPI常量，Fortran 90建议使用MPI模块，但模块未完全支持Fortran 90接口（如带“choice”参数的例程）。  

#### 5. Running Programs with mpiexec  
- **5.1 Standard mpiexec**：介绍MPI标准中的mpiexec参数，如-n指定进程数、-f指定机器文件（格式为主机名:进程数）。  
- **5.2 Extensions for All Process Management Environments**：不同通信子系统和进程管理器的mpiexec扩展参数，力求跨环境统一。  
- **5.3 mpiexec Extensions for the Hydra Process Manager**：Hydra为MPICH默认进程管理器，扩展细节见在线文档。  
- **5.4 Extensions for the gforker Process Management Environment**：gforker用于单节点调试，支持-env设置环境变量、-l标记输出、-maxtime设置超时等参数，也可通过环境变量（如MPIEXEC_TIMEOUT）控制行为。  
- **5.5 Restrictions of the remshell Process Management Environment**：remshell基于ssh启动进程，功能有限，忽略环境变量控制参数，需通过machines文件指定主机。  
- **5.6 Using MPICH with Slurm and PBS**：Hydra原生支持Slurm和PBS，也可通过srun（Slurm）或OSC mpiexec（PBS）启动作业。  

#### 6. Specification of Implementation Details  
- **6.1 MPI Error Handlers for Communicators**：MPICH的MPI通信器错误处理回调函数无需实现特定可变参数，只需包含通信器和错误码参数。  

#### 7. Debugging  
- **7.1 TotalView**：MPICH支持使用TotalView调试器，可通过totalview -a mpiexec启动调试，或使用间接启动功能配置并行参数。  

#### 8. Checkpointing  
- **8.1 Configuring for checkpointing**：需安装BLCR 0.8.2库，通过configure选项--enable-checkpointing --with-hydra-ckpointlib=blcr配置，非默认安装时需指定BLCR路径。  
- **8.2 Taking checkpoints**：使用mpiexec的-ckpointlib指定BLCR，-ckpoint-prefix指定检查点目录，支持手动（SIGUSR1信号）或自动（-ckpoint-interval）触发检查点，可通过-ckpoint-num指定重启的检查点版本。  

#### 9. Other Tools Provided with MPICH  
- **内容**：MPICH包含MPI功能测试套件，位于mpich/test/mpi目录，安装后可通过make testing运行，也可配置用于测试其他MPI实现。  

#### A. Frequently Asked Questions  
- **内容**：常见问题维护于在线文档，链接为https://github.com/pmodels/mpich/blob/main/doc/wiki/faq/Frequently_Asked_Questions.md。  



