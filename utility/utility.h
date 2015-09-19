/*! \file utility.h
    \brief ユーティリティ関数の宣言と実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _UTILITY_H_
#define _UTILITY_H_

#pragma once

#include <boost/cast.hpp>       // for boost::numeric_cast
#include <boost/optional.hpp>   // for boost::optional
#include <boost/utility.hpp>    // for boost::checked_delete

namespace utility {
    template <typename T>
    //! A function.
    /*!
        関数が成功したかどうかを判断する
        \tparam T 関数の戻り値の型
        \param x HRESULTの値
        \return 成功したらboost::optional<HRESULT>、失敗したらboost::none
    */
    boost::optional<HRESULT> v_return(T const & x);

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

    template <typename T> boost::optional<HRESULT> v_return(T const & x)
    {
        auto const hr = boost::numeric_cast<HRESULT>(x);
        return hr >= 0 ? boost::optional<HRESULT>(hr) : boost::none;
    }
}

#endif  // _UTILITY_H_
