/*! \file useavx.h
	\brief AVX命令を使用するかどうかを表す列挙型の宣言

	Copyright ©  2014 @dc1394 All Rights Reserved.
	This software is released under the BSD 2-Clause License.
*/

#ifndef _USEAVX_H_
#define _USEAVX_H_

#pragma once

namespace moleculardynamics {
	enum class UseAVX : bool {
		True = true,
		False = false
	};
}

#endif  // _ELEMENT_H_
