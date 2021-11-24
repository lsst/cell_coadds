// -*- LSST-C++ -*-
/*
 * This file is part of cell_coadds.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef LSST_CELL_COADDS_GridContainer_h
#define LSST_CELL_COADDS_GridContainer_h

#include <iterator>
#include <optional>
#include <vector>

#include "lsst/cell_coadds/GridIndex.h"
#include "lsst/geom/Box.h"

namespace lsst {
namespace cell_coadds {

/**
 *  A semi-private mixin base class for GridContainer and GridContainerBuilder.
 *
 *  _GridContainerCommon injects a few public members into its subclasses, but
 *  otherwise should be considered an implementation detail.  It cannot be
 *  constructed on its own, and is not exposed to Python.
 */
class _GridContainerCommon {
public:
    /// Struct used for indexing.
    using Index = GridIndex;

    /// Return the number of cells in the container in each dimension.
    Index const& get_shape() const { return _shape; }

    /// Return the index of the first cell in the container.
    Index const& get_offset() const { return _offset; }

protected:
    // Construct with shape and default `(0, 0)` offset.
    explicit _GridContainerCommon(Index const& shape) : _shape(shape), _offset{0, 0} {}

    // Construct with shape and explicit offset.
    _GridContainerCommon(Index const& shape, Index const& offset) : _shape(shape), _offset(offset) {}

    // Destructor is protected so it doesn't need to be virtual - nobody should
    // ever construct one of these on its own, or delete through a pointer to
    // one.
    ~_GridContainerCommon() = default;

    // Flatten a 2-d index into a 1-d index.  All derived classes should go
    // through this method to map their 2-d grid to an underlying container.
    std::size_t _flatten_index(Index const& index) const {
        return (index.x - _offset.x) + (index.y - _offset.y) * _shape.x;
    }

    // Return the number of cells in the container.
    std::size_t _size() const { return _shape.x * _shape.y; }

private:
    Index _shape;
    Index _offset;
};

template <typename T>
class GridContainer;

/**
 * A construction helper for GridContainer.
 *
 * @tparam T    Type of objects held by the container.
 *
 * GridContainerBuilder is mutable but is mostly not readable; the expected
 * pattern is to use a temporary builder object briefly to construct a
 * GridContainer that is then used without modification.
 */
template <typename T>
class GridContainerBuilder final : public _GridContainerCommon {
public:
    /// Struct used for indexing.
    using Index = GridIndex;

    /// Construct with shape and default `(0, 0)` offset.
    explicit GridContainerBuilder(Index const& shape) : _GridContainerCommon(shape), _array(_size()) {}

    /// Construct with shape and explicit offset.
    GridContainerBuilder(Index const& shape, Index const& offset)
            : _GridContainerCommon(shape, offset), _array(_size()) {}

    /// Return the number of cells in the container.
    std::size_t size() const { return _array.size(); }

    /**
     * Set the value for a single cell in the container.
     *
     * @param index   2-d index for the cell to set.
     * @param value   Value for the cell.
     *
     * `set` may be called more than once on a cell, but must be called at
     * least once on all cells before `finish` is called.
     */
    void set(Index const& index, T value) { _array[_flatten_index(index)] = std::move(value); }

    /**
     * Finish construction and check that all cells have been set.
     *
     * This consumes the builder in order to move cell values into the new
     * container.  In Python the builder is (shallow) copied instead.
     */
    GridContainer<T> finish() &&;

private:
    // Construct from common geometry with an empty underlying array.
    GridContainerBuilder(_GridContainerCommon const& common)
            : _GridContainerCommon(common), _array(_size()) {}

    template <typename U>
    friend class GridContainer;

    std::vector<std::optional<T>> _array;
};

/**
 * A container whose elements form a 2-d grid.
 *
 * @tparam T    Type of objects held by the container.
 *
 * GridContainer is immutable aside from compiler-generated assignment
 * operators.
 *
 * We expose only an instantiation for `pybind11::object` to Python, allowing
 * it to hold arbitrary Python objects.  We do not provide any access to the
 * assignment operators there, making this a fully immutable container in
 * Python (a more important property there because Python has no `const`).
 * Python type stubs in `_cell_coadds.pyi` allow it to be used as a generic
 * in static typing.
 *
 * The lack of Python wrappers for other instantiations can make it hard to
 * construct C++ objects that use other instantiations internally.  The
 * `rebuild_transformed` method can be used with a callback that uses
 * `pybind11::cast` to work around this (though unfortunately it's difficult
 * to pass an instantiation of `pybind11::cast` directly as the callback due
 * to C++'s templated-function-call parsing rules).
 */
template <typename T>
class GridContainer final : public _GridContainerCommon {
public:
    /// Struct used for indexing.
    using Index = GridIndex;

    /// Immutable iterator type.
    using const_iterator = typename std::vector<T>::const_iterator;

    /// Iterator type.  Always const.
    using iterator = const_iterator;

    /**
     * Construct from a builder instance.
     *
     * @param builder   Builder used to set the values of cells.
     */
    explicit GridContainer(GridContainerBuilder<T> builder) : _GridContainerCommon(builder), _array() {
        _array.reserve(builder.size());
        std::transform(
            std::make_move_iterator(builder._array.begin()),
            std::make_move_iterator(builder._array.end()),
            std::back_inserter(_array),
            [](std::optional<T>&& element) {
                if (!element) {
                    throw LSST_EXCEPT(pex::exceptions::LogicError, "Grid is not complete.");
                }
                return std::move(element.value());
            });
    }

    /**
     * Return an iterator to the first cell.
     *
     * The order of iteration is unspecified, but it is guaranteed to start
     * with one corner of the grid and finish at the opposite corner.
     */
    const_iterator begin() const { return _array.begin(); }

    /// Return an iterator to one past the last cell.
    const_iterator end() const { return _array.end(); }

    /// Return the first cell value.
    T const& get_first() const { return _array.front(); }

    /// Return the last cell value.
    T const& get_last() const { return _array.back(); }

    /// Return the cell value at the given index.
    T const& operator[](Index const& index) const { return _array[_flatten_index(index)]; }

    /// Return the number of cells in the container.
    std::size_t size() const { return _array.size(); }

    /**
     * Create a new `GridContainerBuilder` with the same shape and offset as
     * `this`, with no cell values.
     */
    template <typename U = T>
    GridContainerBuilder<U> rebuild_empty() const {
        return GridContainerBuilder<U>(*this);
    }

    /**
     * Create a new `GridContainerBuilder` with the same shape and offset as
     * `this`, with cell values created by applying a callback function to each
     * cell value in `this`.
     *
     * This consumes the container and moves values into the given function.
     * In Python the container is (shallow) copied instead.
     */
    template <typename F>
    GridContainerBuilder<decltype(std::declval<F>()(std::declval<T>()))> rebuild_transformed(F func) && {
        GridContainerBuilder<decltype(func(std::declval<T>()))> builder(*this);
        std::transform(
            std::make_move_iterator(_array.begin()),
            std::make_move_iterator(_array.end()),
            builder._array.begin(),
            [func](T&& original) { return std::optional(func(std::move(original))); });
        return builder;
    }

private:
    std::vector<T> _array;
};

template <typename T>
GridContainer<T> GridContainerBuilder<T>::finish() && {
    return GridContainer(std::move(*this));
}

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_GridContainer_h
