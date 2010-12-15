/*
 * Copyright 2003-2005 Daniel F. Savarese
 * Copyright 2006-2009 Savarese Software Research Corporation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.savarese.com/software/ApacheLicense-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
 * @file
 * This header defines the kd_tree class and its support classes.
 */

#ifndef __SSRC_SPATIAL_KDTREE_H
#define __SSRC_SPATIAL_KDTREE_H

#include <detail/kd_tree_range_search_iterator.h>
#include <detail/kd_tree_node.h>
#include <rectangle_region.h>

#include <algorithm>
#include <utility>
#include <vector>
#include <limits>

__BEGIN_NS_SSRC_SPATIAL

/**
 * kd_tree_traits stores metadata about kd_tree instances.
 */
template<typename Tree>
struct kd_tree_traits {
  typedef typename Tree::key_type key_type;
  typedef typename Tree::mapped_type mapped_type;
  typedef typename Tree::value_type value_type;
  typedef typename Tree::pointer pointer;
  typedef typename Tree::const_pointer const_pointer;
  typedef typename Tree::reference reference;
  typedef typename Tree::const_reference const_reference;
  typedef typename Tree::discriminator_type discriminator_type;
  typedef typename Tree::node_type node_type;
  typedef typename Tree::iterator iterator;
  typedef typename Tree::const_iterator const_iterator;
  typedef typename Tree::size_type size_type;
  typedef typename key_type::value_type coordinate_type;
  typedef Tree tree_type;

  /**
   * Returns the maximum value a coordinate may assume.
   *
   * @return The maximum value a coordinate may assume.
   */
  static const coordinate_type max_coordinate() {
    return std::numeric_limits<coordinate_type>::infinity();
  }

  /**
   * Returns the minimum value a coordinate may assume.
   *
   * @return The minimum value a coordinate may assume.
   */
  static const coordinate_type min_coordinate() {
    return -std::numeric_limits<coordinate_type>::infinity();
  }    
};


/**
 * kd_tree_const_traits stores metadata about const kd_tree instances.
 */
template<typename Tree>
struct kd_tree_const_traits : public kd_tree_traits<Tree> {
  typedef typename Tree::const_pointer pointer;
  typedef typename Tree::const_reference reference;
};


// Note: we store the discriminator in each node to avoid modulo division,
// trading space for time.
/**
 * A k-d tree divides a k-dimensional space relative to the points it
 * contains by storing them in a binary tree, discriminating by a
 * different dimension at each level of the tree.  It allows efficient
 * point data retrieval (<em>O(lg(n))</em>) and range searching.
 *
 * <p>kd_tree tries to conform to the map interface from the STL where
 * it makes sense.  Therefore, you will find it can be used with some
 * STL algorithms.</p>
 */
template<typename Point,
         typename Value,
         typename Discriminator = unsigned char,
         typename Size = unsigned int>
class kd_tree {
public:

  typedef kd_tree_traits<kd_tree> traits;
  typedef kd_tree_const_traits<kd_tree> const_traits;
  typedef Point key_type;
  typedef Value mapped_type;
  typedef std::pair<const key_type, mapped_type> value_type;
  typedef value_type* pointer;
  typedef pointer const const_pointer;
  typedef value_type& reference;
  typedef const reference const_reference;
  typedef Discriminator discriminator_type;
  typedef rectangle_region<key_type> default_region_type;
  // Is this really what we want--two distinct types as
  // opposed to iterator and const iterator?
  typedef
  detail::kd_tree_range_search_iterator<traits, default_region_type> iterator;
  typedef
  detail::kd_tree_range_search_iterator<const_traits, default_region_type>
  const_iterator;

  typedef Size size_type;

  static const unsigned int dimensions = std::tuple_size<key_type>::value;

  typedef detail::kd_tree_node<traits> node_type;

private:

  struct node_comparator {
    discriminator_type discriminator;

    explicit node_comparator() : discriminator(0) { }

    bool operator()(const node_type* const & n1,
                    const node_type* const & n2) const
    {
      return (n1->point()[discriminator] < n2->point()[discriminator]);
    }
  };

  node_type *_root;
  size_type _size;
  iterator _end_iterator;
  const_iterator _const_end_iterator;

  node_type *
  get_node(const key_type & point, node_type ** const parent = 0) const {
    discriminator_type discriminator;
    node_type *node = _root, *last = 0;

    while(node != 0) {
      discriminator = node->discriminator;

      if(point[discriminator] > node->point()[discriminator]) {
        last = node;
        node = node->child_high;
      } else if(point[discriminator] < node->point()[discriminator]) {
        last = node;
        node = node->child_low;
      } else if(node->point() == point) {
        if(parent != 0)
          *parent = last;
        return node;
      } else {
        last = node;
        node = node->child_high;
      }
    }

    if(parent != 0)
      *parent = last;

    return 0;
  }

  node_type * get_minimum_node(node_type * const node, node_type * const p,
                               const discriminator_type discriminator,
                               node_type ** const parent)
  {
    node_type *result;

    if(discriminator == node->discriminator) {
      if(node->child_low != 0)
        return
          get_minimum_node(node->child_low, node,
                           discriminator, parent);
      else
        result = node;
    } else {
      node_type *nlow = 0, *nhigh = 0;
      node_type *plow, *phigh;

      if(node->child_low != 0)
        nlow = get_minimum_node(node->child_low, node,
                                discriminator, &plow);

      if(node->child_high != 0)
        nhigh = get_minimum_node(node->child_high, node,
                                 discriminator, &phigh);

      if(nlow != 0 && nhigh != 0) {
        if(nlow->point()[discriminator] < nhigh->point()[discriminator]) {
          result  = nlow;
          *parent = plow;
        } else {
          result  = nhigh;
          *parent = phigh;
        }
      } else if(nlow != 0) {
        result  = nlow;
        *parent = plow;
      } else if(nhigh != 0) {
        result  = nhigh;
        *parent = phigh;
      } else
        result  = node;
    }

    if(result == node)
      *parent = p;
    else if(node->point()[discriminator] < result->point()[discriminator]) {
      result  = node;
      *parent = p;
    }

    return result;
  }

  node_type * recursive_remove_node(node_type * const node) {
    discriminator_type discriminator;
    node_type *new_root, *parent;

    if(node->child_low == 0 &&
       node->child_high == 0)
      return 0;
    else
      discriminator = node->discriminator;

    if(node->child_high == 0) {
      node->child_high = node->child_low;
      node->child_low  = 0;
    }

    new_root = get_minimum_node(node->child_high, node,
                                discriminator, &parent);

    if(parent->child_low == new_root)
      parent->child_low = recursive_remove_node(new_root);
    else
      parent->child_high = recursive_remove_node(new_root);

    new_root->child_low  = node->child_low;
    new_root->child_high = node->child_high;
    new_root->discriminator = node->discriminator;

    return new_root;
  }

  // Splitting up remove in this way allows us to implement
  // iterator erase(iterator) properly.
  bool remove(const key_type & point, node_type * const node,
              node_type * const parent) {
    node_type * const new_root = recursive_remove_node(node);

    if(parent == 0)
      _root = new_root;
    else if(node == parent->child_low)
      parent->child_low = new_root;
    else
      parent->child_high = new_root;

    // Must zero children so they are not deleted by ~node_type()
    node->child_low  = 0; 
    node->child_high = 0;

    --_size;
    delete node;

    return true;
  }

  bool add(const key_type & point, const mapped_type & value,
            node_type ** const node, node_type *parent,
            mapped_type * const replaced = 0)
  {
    if(parent == 0) {
      if(_root != 0)
        *node = _root;
      else {
        _root = *node = new node_type(0, point, value);
        ++_size;
        return false;
      }
    } else if(*node == 0) {
      discriminator_type discriminator = parent->discriminator;
      node_type* & child =
        (point[discriminator] >= parent->point()[discriminator] ?
         parent->child_high : parent->child_low);

      if(++discriminator >= dimensions)
        discriminator = 0;

      child = *node = new node_type(discriminator, point, value);

      ++_size;
      return false;
    }

    if(replaced != 0)
      *replaced = (*node)->value();

    (*node)->value() = value;

    return true;
  }
  /*
  template<template<typename> class container>
  static
  node_type * optimize(typename container<node_type*>::iterator begin,
                       typename container<node_type*>::iterator end,
                       node_comparator & comparator)
  */
  template<class container>
  static node_type * optimize(typename container::iterator begin,
                              typename container::iterator end,
                              node_comparator & comparator)
  {
    //typedef typename container<node_type*>::iterator container_iterator;
    typedef typename container::iterator container_iterator;
    node_type *midpoint = 0;
    typename container_iterator::difference_type diff;

    diff = end - begin;

    if(diff > 1) {
      discriminator_type discriminator = comparator.discriminator;
      container_iterator nth = begin + (diff >> 1);
      container_iterator nthprev = nth - 1;

      //std::nth_element(begin, nth, end, comparator);
      std::stable_sort(begin, end, comparator);

      // Ties go in the right subtree.
      while(nth > begin &&
            (*nth)->point()[discriminator] == 
            (*nthprev)->point()[discriminator])
        {
          --nth;
          --nthprev;
        }

      midpoint = *nth;
      midpoint->discriminator = discriminator;

      if(++discriminator >= dimensions)
        discriminator = 0;

      comparator.discriminator = discriminator;

      // Left subtree
      midpoint->child_low =
        optimize<container>(begin, nth, comparator);

      comparator.discriminator = discriminator;

      // Right subtree
      midpoint->child_high =
        optimize<container>(nth + 1, end, comparator);
    } else if(diff == 1) {
      midpoint = *begin;
      midpoint->discriminator = comparator.discriminator;
      midpoint->child_low = 0;
      midpoint->child_high = 0;
    }

    return midpoint;
  }

  /*
  template<template<typename> class container>
  static void fill_container(container<node_type*> & c, node_type *node)
  */
  template<class container>
  static void fill_container(container & c, node_type * const node) {
    if(node == 0)
      return;
    c.push_back(node);
    fill_container(c, node->child_low);
    fill_container(c, node->child_high);
  }

  static void init_point(key_type & point,
                         const typename traits::coordinate_type & value)
  {
    for(unsigned int i=0; i < dimensions; ++i)
      point[i] = value;
  }
  
  static const key_type upper_bound() {
    key_type bound;
    init_point(bound, traits::max_coordinate());
    return bound;
  }

  static const key_type lower_bound() {
    key_type bound;
    init_point(bound, traits::min_coordinate());
    return bound;
  }

public:

  // TODO: need assignment operator.

  /**
   * Constructs an empty kd_tree.
   */
  explicit kd_tree() :
    _root(0), _size(0), _end_iterator(), _const_end_iterator()
  { }

  /**
   * Constructs a copy of a kd_tree.
   *
   * @param tree The kd_tree to copy.
   */
  kd_tree(const kd_tree & tree) :
    _root(0), _size(0), _end_iterator(), _const_end_iterator()
  {
    for(const_iterator p = tree.begin(); !p.end_of_range(); ++p)
      insert(p->first, p->second);
  }

  /**
   * Deallocates the memory allocated by the kd_tree.
   */
  virtual ~kd_tree() { delete _root; }

  /**
   * Removes all elements from the container, leaving it empty.
   */
  void clear() {
    delete _root;
    _root = 0;
    _size = 0;
  }

  /**
   * Returns the number of point-value mappings in the kd_tree.
   *
   * @return The number of point-value mappings in the kd_tree.
   */
  const size_type size() const {
    return _size;
  }

  /**
   * Returns the maximum number of elements that can be contained.
   *
   * @return The maximum number of elements that can be contained.
   */
  const size_type max_size() const {
    return std::numeric_limits<size_type>::max();
  }

  /**
   * Returns true if the container has no elements, false if it
   * contains one or more elements.
   *
   * @return true if the container has no elements, false if it
   * contains one or more elements.
   */
  bool empty() const {
    return (_root == 0);
  }

  /**
   * Returns an iterator pointing to the leftmost bottommost point in
   * the container.  If empty, it equals end().
   *
   * @return An iterator pointing to the leftmost bottommost point in
   * the container.  If empty, it equals end().
   */
  iterator begin() {
    return iterator(default_region_type(lower_bound(), upper_bound()), _root);
  }

  /**
   * Returns a const iterator pointing to the leftmost bottommost point in
   * the container.  If empty, it equals end().
   *
   * @return A const iterator pointing to the leftmost bottommost point in
   * the container.  If empty, it equals end().
   */
  const_iterator begin() const {
    return const_iterator(default_region_type(lower_bound(), upper_bound()),
                          _root);
  }

  /**
   * Returns an iterator denoting the end of a range.
   *
   * @return An iterator denoting the end of a range.
   */
  iterator & end() {
    return _end_iterator;
  }

  /**
   * Returns a const iterator denoting the end of a range.
   *
   * @return A const iterator denoting the end of a range.
   */
  const const_iterator & end() const {
    return _const_end_iterator;
  }

  /**
   * Inserts a key value pair, replacing any existing value that may
   * be present, and optionally retrieving the replaced value.
   *
   * @param point The key corresponding to the value to be inserted
   * @param value The value to be inserted.
   * @param replaced If you want to retrieve a value that was
   * replaced, have this parameter point to a valid object to store
   * the value.  If not, its default value of zero will prevent it
   * from being retrieved.
   * @return true if value is replaced, false if not.
   */
  bool insert(const key_type & point, const mapped_type & value,
              mapped_type * const replaced = 0)
  {
    node_type *parent;
    node_type *node = get_node(point, &parent);

    return add(point, value, &node, parent, replaced);
  }

  /**
   * Inserts a key value pair, but--in conformance with the STL
   * container interface--doesn't replace any existing value that may
   * be present.  Note, this method is slower than
   * bool insert(const key_type &, const mapped_type, mapped_type*)
   *
   * @return A pair containing the iterator pointing to the inserted
   * value and true if the mapping is inserted.  If the value is not
   * inserted, a pair containing the iterator pointing to the existing
   * value and false.
   */
  std::pair<iterator,bool> insert(const value_type & mapping) {
    // Ideally, we'd do this all in one step, but that will have
    // to wait until we optimize the way we handle iterators.
    mapped_type existing;
    const bool replaced  = insert(mapping.first, mapping.second, &existing);
    const iterator value = find(mapping.first);

    if(replaced)
      value._node->value() = existing;

    return std::pair<iterator,bool>(value,!replaced);
  }


  /**
   * Retrieves the value at the given location.  In conformance with
   * the STL map interface, if the point key does not occur in the
   * kd_tree, a new value is inserted and returned using the default
   * constructor for kd_tree::mapped_type.
   *
   * @param point The location from which to retrieve the value.
   * @return The value at the given location.
   */
  mapped_type & operator[](const key_type & point) {
    node_type *parent;
    node_type *node = get_node(point, &parent);

    if(node == 0)
      add(point, mapped_type(), &node, parent);

    return node->value();
  }

  /**
   * Removes the point-value mapping corresponding to the given point key.
   * Optionally, retrieves the removed value if a mapping was present.
   *
   * @param point The point key of the mapping to remove.
   * @param erased A pointer to a mapped_type in which to store the
   * erased value.  This pointer may be null if you do not want to
   * retrieve the removed value and only want to remove the mapping.
   * @return true if a mapping was present and removed, false if not.
   */
  bool remove(const key_type & point, mapped_type * const erased = 0) {
    node_type *parent;
    node_type * const node = get_node(point, &parent);

    if(node == 0)
      return false;

    if(erased != 0)
      *erased = node->value();

    return remove(point, node, parent);
  }

  /**
   * Removes the point-value mapping corresponding to the given point key.
   *
   * @param point The point key of the mapping to remove.
   * @return The number of mappings erased (0 or 1).
   */
  size_type erase(const key_type & point) {
    return remove(point);
  }

  /**
   * Removes the point-value mapping at the specified position.  The
   * supplied iterator becomes invalid upon a successful erasure.
   *
   * @param pos The kd_tree::iterator denoting the location of the mapping.
   * @return An iterator pointing to the next element after pos in the
   * range search sequence containing pos.  If no such element exists,
   * or if the element pointed to by pos is not contained in the tree
   * (making it impossible to remove from the tree), end() is returned.
   */
  iterator erase(iterator pos) {
    if(pos.end_of_range())
      return _end_iterator;

    node_type *parent;
    node_type * const node = get_node(pos->first, &parent);
    
    if(node == 0)
      return _end_iterator;

    typename iterator::stack_type & stack = pos._stack;

    // Pop any children.  Tree at parent and above is unchanged.
    // Low child is pushed last so check it first.
    if(!stack.empty() && node->child_low == stack.top()) {
      stack.pop();
    }
    if(!stack.empty() && node->child_high == stack.top()) {
      stack.pop();
    }

    const bool low_child = (parent != 0 && parent->child_low == node);

    if(remove(pos->first, node, parent)) {
      if(parent != 0) {
        if(low_child && parent->child_low != 0) {
          stack.push(parent->child_low);
        } else if(!low_child && parent->child_high != 0) {
          stack.push(parent->child_high);
        }
        pos.advance();
        return pos;
      } else if(_root != 0) {
        stack.push(_root);
        pos.advance();
        return pos;
      }
    }

    return _end_iterator;
  }

  /**
   * Returns an iterator for mappings that are contained in the
   * rectangle defined by the given lower left-hand and upper
   * right-hand corners.  The mappings returned include those occuring
   * at points on the bounding rectangle, not just those inside.
   *
   * @param lower The lower left-hand corner of the bounding
   * rectangle.
   * @param upper The upper right-hand corner of the bounding
   * rectangle.
   * @return A kd_tree::iterator for mappings that are contained in the
   * specified rectangle.
   */
  iterator begin(const key_type & lower, const key_type & upper) {
    return iterator(default_region_type(lower, upper), _root);
  }

  /**
   * Returns a const iterator for mappings that are contained in the
   * rectangle defined by the given lower left-hand and upper
   * right-hand corners.  The mappings returned include those occuring
   * at points on the bounding rectangle, not just those inside.
   *
   * @param lower The lower left-hand corner of the bounding
   * rectangle.
   * @param upper The upper right-hand corner of the bounding
   * rectangle.
   * @return A kd_tree::const_iterator for mappings that are contained in the
   * specified rectangle.
   */
  const_iterator begin(const key_type & lower, const key_type & upper) const {
    return const_iterator(default_region_type(lower, upper), _root);
  }

  // TODO: Document these and rectangle_region.  Implement circle_region
  // and sphere_region and write unit tests.  Implement nearest neighbor
  // search and k-nearest neighbors search.  Try to implement nearest
  // neighbor as standalone generic functions.
  // Move kd_tree_range_search_iterator out of detail.
  template<typename Region>
  detail::kd_tree_range_search_iterator<traits, Region>
  begin(const Region & region) {
    return
      detail::kd_tree_range_search_iterator<traits, Region>(region, _root);
  }

  template<typename Region>
  detail::kd_tree_range_search_iterator<traits, Region> end() {
    return detail::kd_tree_range_search_iterator<traits, Region>();
  }

  template<typename Region>
  detail::kd_tree_range_search_iterator<const_traits, Region>
  begin(const Region & region) const {
    return
      detail::kd_tree_range_search_iterator<const_traits, Region>(region, _root);
  }

  template<typename Region>
  detail::kd_tree_range_search_iterator<const_traits, Region> end() const {
    return detail::kd_tree_range_search_iterator<const_traits, Region>();
  }

  /**
   * Returns true if the kd_tree contains a point-value mapping for the
   * specified point, false if not.  Optionally, retrieves the value
   * at the given location (faster than find()).
   *
   * @param point The location from which to retrieve the value.
   * @param value A pointer to a kd_tree::mapped_type in which to store the
   * retrieved value.  This pointer may be null if you do not want to
   * retrieve the value and only want to see if the point is present in
   * the kd_tree.
   * @return true if the kd_tree contains a point-value mapping for the
   * specified point, false if not.
   */
  bool get(const key_type & point, mapped_type * const value = 0) const {
    const node_type * const node = get_node(point);

    if(node == 0)
      return false;
    else if(value != 0)
      *value = node->value();

    return true;
  }

  /**
   * Returns the location of the point-value mapping for the specified
   * point (slower than get()).  If there is no mapping present, the
   * iterator is equivalent to end().
   *
   * @param point The key of the poinit-value mapping to find.
   * @return true A kd_tree::iterator for the location of the
   * point-value mapping matching the specified point.  end() if there
   * is no mapping.
   */
  iterator find(const key_type & point) {
    return iterator(default_region_type(point, upper_bound()), _root, true);
  }

  /**
   * Returns the location of the point-value mapping for the specified
   * point (slower than get()).  If there is no mapping present, the
   * iterator is equivalent to end().
   *
   * @param point The key of the poinit-value mapping to find.
   * @return true A kd_tree::const_iterator for the location of the
   * point-value mapping matching the specified point.  end() if there
   * is no mapping.
   */
  const_iterator find(const key_type & point) const {
    return const_iterator(default_region_type(point, upper_bound()),
                          _root, true);
  }

  // Balances the tree.  Very expensive!
  /**
   * Optimizes the performance of future search operations by balancing the
   * kd_tree.  The balancing operation is relatively expensive, but can
   * significantly improve the performance of searches.  Usually, you
   * don't have to optimize a tree which contains random key values
   * inserted in a random order.
   */
  void optimize() {
    if(empty())
      return;

    typedef std::vector<node_type*> container;
    container nodes;

    nodes.reserve(size());
    //fill_container<std::vector>(nodes, _root);
    fill_container(nodes, _root);

    node_comparator comparator;
    //_root = optimize<std::vector>(nodes.begin(), nodes.end(), comparator);
    _root = optimize<container>(nodes.begin(), nodes.end(), comparator);
  }

};


__END_NS_SSRC_SPATIAL

#endif
