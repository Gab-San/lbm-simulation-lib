#pragma once

#include "lbm-sim/core/point.hpp"
#include "lbm-sim/core/vector.hpp"

#include "quill/DeferredFormatCodec.h"
#include "quill/bundled/fmt/ostream.h"

#include "quill/core/Codec.h"

template <typename T, unsigned short int dim>
struct fmtquill::formatter<lbm::utils::Point<T, dim>>
    : fmtquill::ostream_formatter {};

template <typename T, unsigned short int dim>
struct quill::Codec<lbm::utils::Point<T, dim>>
    : quill::DeferredFormatCodec<lbm::utils::Point<T, dim>> {};

template <typename T, unsigned short int dim>
struct fmtquill::formatter<lbm::utils::Vector<T, dim>>
    : fmtquill::ostream_formatter {};

template <typename T, unsigned short int dim>
struct quill::Codec<lbm::utils::Vector<T, dim>>
    : quill::DeferredFormatCodec<lbm::utils::Vector<T, dim>> {};
