/*! \file utility.h
    \brief ユーティリティ関数の宣言と実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _UTILITY_H_
#define _UTILITY_H_

#pragma once

#include <stdexcept>            // for std::runtime_error
#include <string>               // for std::to_string
#include <boost/cast.hpp>       // for boost::numeric_cast
#include <boost/utility.hpp>    // for boost::checked_delete

namespace utility {
    template <typename T>
    //! A function.
    /*!
        関数が成功したかどうかを判断する
        \tparam T 関数の戻り値の型
        \param x HRESULTの値
    */
    void v_return(T const & x);

    template <typename T>
    //! A struct.
    /*!
        リソースを安全に解放するクラス
        \tparam T リソースの型
    */
    struct Safe_Release {
        //! A public member function.
        /*!
            リソースを安全に解放する
            \param p リソースへのポインタ
        */
        void operator()(T * p) {
            if (p) {
                p->Release();
                p = nullptr;
            }
        }
    };

    template <typename T>
    //! A struct.
    /*!
        確保したメモリを安全に解放するクラス
        \tparam T 確保したメモリの型
    */
    struct Safe_Delete {
        //! A public member function.
        /*!
            確保したメモリを安全に解放する
            \param p 確保したメモリの先頭アドレス
        */
        void operator()(T * p) {
            if (p) {
                boost::checked_delete(p);
                p = nullptr;
            }
        }
    };

    template <typename T> void v_return(T const & x)
    {
        auto const hr = boost::numeric_cast<HRESULT>(x);
        if (hr < 0) {
            throw std::runtime_error("function Failed! HRESULT: " + std::to_string(hr));
        }
    }
}

#endif  // _UTILITY_H_
