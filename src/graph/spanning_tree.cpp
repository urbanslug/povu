#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./spanning_tree.hpp"
#include "bracket_list.hpp"


namespace povu::spanning_tree {
using namespace povu::constants;
using namespace povu::graph_types;



/*
 * Edge
 * ----
 */
// constructor(s)
Edge::Edge(): null_(true) {}

Edge::Edge(std::size_t id, std::size_t src, std::size_t tgt, color c):
  id_(id), src(src), tgt(tgt), null_(false), color_(c) {}

// getters
std::size_t Edge::id() const { return this->id_; }
color Edge::get_color() const { return this->color_; }
std::size_t Edge::get_parent() const { return this->src; }
std::size_t Edge::get_parent_v_idx() const { return this->src; }
std::size_t Edge::get_v1() const { return this->src; }
std::size_t Edge::get_v2() const { return this->tgt; }
std::size_t Edge::get_child() const { return this->tgt; }
std::size_t Edge::get_class() const { return this->class_; }
std::size_t Edge::get_class_idx() { return this->class_; }

// setters
void Edge::set_class_idx(std::size_t c) { this->class_ = c; }
void Edge::set_class(std::size_t c) { this->class_ = c; }


/*
 * BackEdge
 * --------
 */

BackEdge::BackEdge(std::size_t id,  std::size_t src,  std::size_t tgt, EdgeType t, pgt::color c)
  : id_(id), src(src), tgt(tgt), class_(INVALID_ID), recent_class_(INVALID_ID),
    recent_size_(INVALID_ID), type_(t), null_(false), color_(c) {}

std::size_t BackEdge::id() const { return this->id_; }
std::size_t BackEdge::get_src() const { return this->src; }
std::size_t BackEdge::get_tgt() const { return this->tgt; }

std::size_t BackEdge::get_class() const { return this->class_; }

color BackEdge::get_color() const { return this->color_; }

bool BackEdge::is_class_defined() const { return this->class_ != INVALID_ID; }
bool BackEdge::is_capping_backedge() const { return this->type_ == EdgeType::capping_back_edge; }

void BackEdge::set_class(std::size_t c) { this->class_ = c; }
EdgeType BackEdge::type() const { return this->type_; }


/*
 * Vertex
 * ------
 */
// Constructor(s)

Vertex::Vertex(std::size_t dfs_num, const std::string &name, VertexType type_)
  : dfs_num_(dfs_num), parent_id(INVALID_ID), name_(name), type_(type_),
    hi_(std::numeric_limits<size_t>::max()), null_(false){}


// getters
std::string const& Vertex::name() const { return this->name_; }
VertexType Vertex::type() const { return this->type_; }
std::size_t  Vertex::hi() const { return this->hi_; }
std::size_t  Vertex::dfs_num() const { return this->dfs_num_; }
bool Vertex::is_leaf() const { return this->children.empty(); }
std::size_t  Vertex::parent() const { return this->parent_id; }
std::set<size_t> const &Vertex::get_ibe() const { return this->ibe; }
std::set<size_t> const &Vertex::get_obe() const { return this->obe; }
size_t const& Vertex::get_parent_idx() const { return this->parent_id; }
size_t Vertex::get_parent_e_idx() const { return this->parent_id; }
std::set<size_t> const& Vertex::get_children() const { return this->children; }
bool Vertex::is_root() const {
  return !this->null_ && this->parent_id == INVALID_IDX;
}

bool Vertex::is_null() const { return this->null_; }

// setters
void Vertex::unset_null(){ this->null_ = false;}
void Vertex::add_obe(std::size_t obe_id) { this->obe.insert(obe_id); }
void Vertex::add_ibe(std::size_t ibe_id) { this->ibe.insert(ibe_id); }
void Vertex::add_child(std::size_t e_id) { this->children.insert(e_id); }
void Vertex::set_parent(std::size_t n_id) { this->parent_id = n_id; }
void Vertex::set_name(std::string const& name) { this->name_ = name; }
void Vertex::set_type(VertexType t) { this->type_ = t; }
void Vertex::set_hi(std::size_t val) { this->hi_ = val; }
void Vertex::set_dfs_num(std::size_t idx) { this->dfs_num_ = idx; }

/*
 * Tree
 * ----
 */

// Constructor(s)

Tree::Tree(std::size_t size) :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  //bracket_lists(std::vector<BracketList*>{}),
  sort_(std::vector<std::size_t>{}),
  sort_g(std::vector<std::size_t>{}),
  equiv_class_count_(0) {
  this->nodes.reserve(size);
  this->tree_edges.reserve(size);
  this->back_edges.reserve(size);
  // this->bracket_lists.reserve(size);
  // TODO: is this necessary?
  this->bracket_lists = std::vector<WBracketList*>(size, nullptr);
}

Tree::~Tree(){
  this->nodes.clear();
  this->tree_edges.clear();
  this->back_edges.clear();
  for (std::size_t i = 0; i < this->bracket_lists.size(); ++i) {
    if (this->bracket_lists.at(i) != nullptr) {
      delete this->bracket_lists[i];
    }
    //else {
    //  std::cout << "Bracket list " << i << " is null" << std::endl;
    //}
  }
  this->bracket_lists.clear();

  this->sort_.clear();
  this->sort_g.clear();
}

void Tree::set_sort(std::size_t idx, std::size_t vertex) {
  this->sort_.at(idx) = vertex;
}

void Tree::set_sort_g(std::size_t idx, std::size_t vertex) {
  this->sort_g.at(idx) = vertex;
}

void Tree::set_dfs_num(std::size_t vertex, std::size_t dfs_num) {
  this->nodes.at(vertex).set_dfs_num(dfs_num);
}

void Tree::set_vertex_type(std::size_t vertex, VertexType type) {
  this->nodes.at(vertex).set_type(type);
}

void Tree::add_vertex(Vertex&& v) {
  this->nodes.push_back(v);
}

Vertex& Tree::get_root()  { return this->nodes.at(this->get_root_idx()); }
std::size_t Tree::get_root_idx() const { return this->root_node_index; }


std::size_t Tree::size() const { return this->nodes.size(); }
std::size_t Tree::tree_edge_count() const { return this->tree_edges.size(); }
std::size_t Tree::back_edge_count() const { return this->back_edges.size(); }

Vertex const &Tree::get_vertex(std::size_t vertex) const {
  return this->nodes.at(vertex);
}

Vertex& Tree::get_vertex_mut(std::size_t vertex) {
  return this->nodes.at(vertex);
}


Vertex const& Tree::get_p_vtx(std::size_t v_idx) const {
  std::size_t p_idx = this->get_parent_v_idx(v_idx);
  return this->get_vertex(p_idx);
}


std::size_t Tree::list_size(std::size_t vertex) {
  return this->bracket_lists.at(vertex)->size();
}

std::size_t Tree::get_hi(std::size_t vertex) {
  return this->nodes.at(vertex).hi();
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_children_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    res.insert(
      std::make_pair(
        this->tree_edges.at(e_idx).id(),
        this->tree_edges.at(e_idx).get_child()
      )
   );
  }
  return res;
}

std::vector<Edge> Tree::get_child_edges(std::size_t vertex) {
  std::vector<Edge> v{};

  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    v.push_back(this->tree_edges.at(e_idx));
  }

  return v;
}

Edge const& Tree::get_parent_edge(std::size_t vertex) const {
  return this->tree_edges.at(this->nodes.at(vertex).get_parent_idx());
}

std::set<std::size_t> Tree::get_obe_idxs(std::size_t vertex) {
  std::set<std::size_t> v {};

  for (auto be_idx : this->nodes.at(vertex).get_obe()) {
    v.insert(be_idx);
  }

  return v;
}

std::set<std::size_t> Tree::get_ibe_idxs(std::size_t vertex) {
  std::set<std::size_t> v{};

  for (auto be_idx : this->nodes.at(vertex).get_ibe()) {
    v.insert(be_idx);
  }

  return v;
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_obe_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_obe()) {
    res.insert(
      std::make_pair(
        this->back_edges.at(e_idx).id(),
        this->back_edges.at(e_idx).get_tgt())
    );
  }
  return res;
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_ibe_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
    res.insert(
      std::make_pair(
        this->back_edges.at(e_idx).id(),
        this->back_edges.at(e_idx).get_src())
    );
  }
  return res;
}

std::set<std::size_t> Tree::get_children(std::size_t vertex) {
  std::set<std::size_t> res {};
  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    res.insert(this->tree_edges.at(e_idx).get_child());
  }
  return res;
}

std::set<std::size_t> Tree::get_ibe(std::size_t vertex) {
  std::set<std::size_t> res {};
  for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
    res.insert(this->back_edges.at(e_idx).get_src());
  }
  return res;
}

std::set<std::size_t> Tree::get_obe(std::size_t vertex) {
  std::set<std::size_t> res{};
  for (auto e_idx : this->nodes.at(vertex).get_obe()) {
    res.insert(this->back_edges.at(e_idx).get_tgt());
  }
  return res;
}

bool Tree::is_root(std::size_t vertex) const {
  return this->get_vertex(vertex).is_root();
}

bool Tree::is_leaf(std::size_t vertex) const {
  return this->get_vertex(vertex).is_leaf();
}

bool Tree::has_child(std::size_t vertex, std::size_t child_idx)  {
  return this->get_children(vertex).count(child_idx);
}

bool Tree::has_ibe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_ibe(vertex).count(qry_idx);
}

bool Tree::has_obe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_obe(vertex).count(qry_idx);
}

Edge& Tree::get_incoming_edge(std::size_t vertex) {
  std::size_t e_idx = this->nodes.at(vertex).get_parent_idx();
  return this->tree_edges.at(e_idx);
}

std::size_t Tree::get_parent(std::size_t vertex) {
  std::size_t i = this->get_vertex(vertex).get_parent_idx();
  return this->tree_edges.at(i).get_parent();
}

std::size_t Tree::get_parent_v_idx(std::size_t v_idx) const {
  std::size_t e_idx = this->get_vertex(v_idx).get_parent_e_idx();
  return this->tree_edges.at(e_idx).get_parent();
}

const Edge& Tree::get_tree_edge(std::size_t edge_idx) const {
    return this->tree_edges.at(edge_idx);
}

std::size_t Tree::get_graph_edge_id(std::size_t tree_edge_id) const {
    return this->tree_graph_idx_map_.at(tree_edge_id);
}

const std::pair<EdgeType, std::size_t>& Tree::get_edge_idx(std::size_t edge_id) const {
  return this->edge_id_map_.at(edge_id);
}

BackEdge& Tree::get_backedge(std::size_t backedge_idx) {
  return this->back_edges.at(backedge_idx);
}

BackEdge &Tree::get_backedge_ref_given_id(std::size_t backedge_id) {
  std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
  return this->back_edges.at(be_idx);
}

BackEdge Tree::get_backedge_given_id(std::size_t backedge_id) {
  std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
  return this->back_edges[be_idx];
}

void Tree::add_tree_edge(std::size_t frm, std::size_t to, std::size_t g_edge_idx, color c) {
  std::size_t edge_idx = this->tree_edges.size();
  std::size_t edge_count = edge_idx + this->back_edges.size();
  this->tree_edges.push_back(Edge(edge_count, frm, to, c));

  //if (this->g_edge_idx_map.count(g_edge_idx)) {
    //std::cerr << "TE g_e id" << g_edge_idx << " id " << edge_count << "\n";
        //throw std::runtime_error("Tree Edge already exists");
  //}

  this->g_edge_idx_map[g_edge_idx] = std::make_pair(EdgeType::tree_edge, edge_idx);

  this->tree_graph_idx_map_[edge_count] = g_edge_idx;

  this->edge_id_map_[edge_count] = std::make_pair(EdgeType::tree_edge, edge_idx);

  this->nodes[frm].unset_null();
  this->nodes[to].unset_null();

  this->nodes[frm].add_child(edge_idx);
  this->nodes[to].set_parent(edge_idx);
}

std::size_t Tree::add_be(std::size_t frm, std::size_t to, EdgeType t, color c) {
  std::size_t back_edge_idx = this->back_edges.size();
  std::size_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(BackEdge(edge_count, frm, to, t, c));
  this->nodes[frm].add_obe(back_edge_idx);
  this->nodes[to].add_ibe(back_edge_idx);

  this->be_id_to_idx_map_[edge_count] = back_edge_idx;

  return back_edge_idx;
}

std::size_t Tree::add_be(std::size_t frm, std::size_t to, std::size_t g_edge_id, EdgeType t, color c) {
  std::size_t back_edge_idx = this->back_edges.size();
  std::size_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(BackEdge(edge_count, frm, to, t, c));

  if (t != EdgeType::capping_back_edge) {
    this->tree_graph_idx_map_[edge_count] = g_edge_id;
    this->g_edge_idx_map[g_edge_id] = std::make_pair(EdgeType::back_edge, back_edge_idx);
  }

  this->edge_id_map_[edge_count] = std::make_pair(EdgeType::back_edge, back_edge_idx);

  this->nodes[frm].add_obe(back_edge_idx);
  this->nodes.at(to).add_ibe(back_edge_idx);


  this->be_id_to_idx_map_[edge_count] = back_edge_idx;

  return back_edge_idx;
}


void Tree::set_hi(std::size_t vertex, std::size_t val) {
  this->nodes.at(vertex).set_hi(val);
}

/**
 * insert the elements of the child bracket list at the
 * beginning of the parent bracket list
 *
 * TODO: make constant time
 *
 * @param parent_vertex
 * @param child_vertex
*/
void Tree::concat_bracket_lists(std::size_t parent_vertex, std::size_t child_vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::Tree::{}]", __func__);

  WBracketList* bl_p = this->bracket_lists[parent_vertex];
  WBracketList* bl_c = this->bracket_lists[child_vertex];

  if (bl_p == nullptr) {
    this->bracket_lists[parent_vertex] = bl_c;
    this->bracket_lists[child_vertex] = nullptr;
  }
  else {
    bl_p->concat(bl_c);
  }
}

// TODO: once deleted do we care to reflect changes in the concated ones?
// WE DO!!
/**
 * delete the bracket that is associated with the backedge
 * given a vertex id and a backedge idx
 */
void Tree::del_bracket(std::size_t vertex, std::size_t backedge_idx) {
  std::string fn_name = std::format("[povu::spanning_tree::Tree::{}]", __func__);

  std::size_t be_id = this->back_edges.at(backedge_idx).id();
  this->bracket_lists[vertex]->del(be_id);
}


void Tree::push(std::size_t vertex, std::size_t backege_idx) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);

  // TODO: based on the Tree constructor we expect the pointer at v_idx will
  // never be null why then do we need to check for null else code fails
  // we then create a  bracket using the backedge ID

  if (this->bracket_lists[vertex] == nullptr) {
    this->bracket_lists[vertex] = new WBracketList{};
  }

  this->bracket_lists[vertex]->push(Bracket(this->back_edges.at(backege_idx).id()));
}


BracketList& Tree::get_bracket_list(std::size_t vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);
  if (this->bracket_lists[vertex] == nullptr) {
    throw std::runtime_error(std::format("{} Bracket list is null", fn_name));
  }

  return this->bracket_lists[vertex]->get_bracket_list();
}


Bracket& Tree::top(std::size_t vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);

  return this->bracket_lists[vertex]->top();
}


std::size_t Tree::new_class() { return this->equiv_class_count_++; }


/**
 * Get the node id of the node with the given sort value
 *
 * The sort vector is sorted in ascending order based on the index from 0..n
 * the value at index zero in the sort vector (`sort_`) is the node id of the
 * node with the smallest sort value and so forth.
 * We can then use this node id to get the node from the nodes vector (`nodes`)
 *
 * @param[in] idx the sort value of the node
 * @return the node id of the node with the given sort value
 */
std::size_t Tree::get_sorted(std::size_t idx) { return  this->sort_.at(idx);}

std::size_t Tree::get_sorted_g(std::size_t idx) { return  this->sort_g.at(idx);}

const std::map<std::size_t, std::pair<EdgeType, std::size_t>>& Tree::get_g_edge_idx_map() const {
  return this->g_edge_idx_map;
}

void Tree::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = LR;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (std::size_t i{}; i < this->size(); i++){
    std::string v_type_str = this->get_vertex(i).type() == VertexType::l ? "-" : "+";
    std::cout << std::format("\t{} [label =  \"{} ({}{})\"];\n", i, i, this->get_vertex(i).name(), v_type_str);
  }

  for (std::size_t i{}; i < this->size(); i++) {
    for (auto c : this->get_child_edges(i)) {
      std::string cl = c.get_class_idx() > 10000
                           ? "\u2205"
                           : std::to_string(c.get_class_idx());

      std::string color = c.get_color() == color::gray ? "gray" : "black";

      std::cout << std::format("\t{}  -- {}  [label=\"{} {}\" color=\"{}\"];\n",
                               i, c.get_child(), c.id(), cl, color);
    }

    for (auto o: this->get_obe_w_id(i)) {
      std::string cl = o.first > 10000 ? "\u2205" : std::to_string(o.first);
      // a capping backedge is red and can never have been gray
      BackEdge be = this->get_backedge_given_id(o.first);

      bool is_capping = be.is_capping_backedge();

      std::string class_ = be.get_class() == UNDEFINED_SIZE_T ?  "" : std::to_string(be.get_class());

      std::string color{};

      if (is_capping) {
        color = "red";
      } else if (this->get_backedge_given_id(o.first).get_color() == color::gray) {
        color = "gray";
      } else if (this->get_backedge_given_id(o.first).get_color() == color::black) {
        color = "black";
      } else {
        color = "blue";
      }

      std::cout << std::format(
          "\t{} -- {} [label=\"{} {}\" style=\"dotted\" penwidth=\"3\" color=\"{}\"];\n",
          i, be.get_tgt(), cl, class_, color
          );
    }
  }

  std::cout << "}" << std::endl;
}

} // namespace spanning_tree
