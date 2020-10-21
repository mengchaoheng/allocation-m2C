1.matlab的安装路径中的toolbox文件夹已经包含QCAT文件夹，并且matlab的默认路径包含该文件及其子文件夹，其余地方的函数仅供参考或修改后用。
2.将qcat1_2_1的内核函数修改成各种不同用途的m代码，再进一步转c代码。注意内存问题，内存不足会导致程序卡死。
3. 生成示例以及修改说明：/*
 * File: dir_alloc_mch.c
 * 参考Efficient Active Set Algorithms for Solving Constrained Least Squares Problems in Aircraft Control Allocation
 * 参考网址 http://research.harkegard.se/qcat/
 * 由matlab代码生成
	Generated on:	23-Nov-2019 12:19:25
	Build type:	Static Library
	Output file:	C:\Users\mengc\Desktop\????matlab?C??_11_20\????\codegen\lib\wls_alloc_mch\wls_alloc_mch.lib
	Processor:	ARM Compatible->ARM Cortex
	Version:	MATLAB Coder 4.0 (R2018a)
 * 生成的代码单独在一个c文件——控制分配实现文件（如wls_alloc_m.c），控制分配实现函数（如wls_alloc_m）将在allocation.c文件——控制分配文件
 * 的control_allocation函数——控制分配函数进行调用
 * allocation.c文件包含control_allocation_initialize用以初始化控制分配结构体，控制分配函数control_allocation
 * 作为控制分配的主调函数。重新生成只需要在allocation.h修改/添加控制分配实现头文件（如wls_alloc_m.h）
 * 最后在调需要的地方调用 control_allocation_initialize，和 control_allocation
 */