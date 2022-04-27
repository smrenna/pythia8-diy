#ifndef DIYGEN_Reduce_HPP
#define DIYGEN_Reduce_HPP 1

// --- callback functions ---//
//
// callback function for merge operator, called in each round of the reduction
// one block is the root of the group
// link is the neighborhood of blocks in the group
// root block of the group receives data from other blocks in the group and reduces the data
// nonroot blocks send data to the root
//

// have not yet done a local neighbor exchange yet.
// this will become a template with parameter that is the Block type

#include "config.hpp"

template <typename B>
void reduceData(B* b,
	    const diy::ReduceProxy& rp,
	    const diy::RegularMergePartners& partners)
{
  typedef decltype(std::declval<B>().reduceData()) return_type;
  typedef typename std::remove_reference<return_type>::type const_type;
  typedef typename std::remove_const<const_type>::type data_type;
  //unsigned round = rp.round(); // current round number
  // can accumulate into the local variable and then enqueue that.
  if(b->reduceBuffer().size()!=b->reduceData().size())
    b->reduceBuffer() = b->reduceData();
  
  data_type& in_vals(b->reduceBuffer()); //reduceData());

  // step 1: dequeue and merge
  for (int i=0; i < rp.in_link().size(); ++i)
    {
      int nbr_gid = rp.in_link().target(i).gid;      
      if (nbr_gid == rp.gid()) continue;

      data_type tmp;
      rp.dequeue(nbr_gid, tmp);
      typename data_type::value_type h = in_vals[0];
      B::reduce(in_vals, tmp);
      //fmt::print("dequeue action {}: {} {} {} {} {}\n",rp.gid(),nbr_gid,h,tmp[0], in_vals[0], round);
    }

  // step 2: enqueue
  for (int i=0; i < rp.out_link().size(); ++i) // redundant since size should equal to 1
    {
      // only send to root of group, but not self
      if (rp.out_link().target(i).gid != rp.gid())
	{
	  rp.enqueue(rp.out_link().target(i), in_vals); // b->reduceData());
	  //fmt::print("equeue action {}: {} {}\n",rp.gid(),in_vals[0], round);
	}
    }
  //b->reduceBuffer().swap(in_vals);
}


template <typename B>
void bc_pointconfig(B* b,                                  // local block
         const diy::ReduceProxy& rp,                // communication proxy
         const diy::RegularBroadcastPartners& partners) // partners of the current block
{
    //unsigned   round    = rp.round();               // current round number
    // step 1: dequeue
    for (int i=0; i < rp.in_link().size(); ++i)
    {
      int nbr_gid = rp.in_link().target(i).gid;
      if (nbr_gid == rp.gid()) continue;

      PointConfig  tmp;
      rp.dequeue(nbr_gid, tmp);
      b->state=tmp;
    }

    if (rp.out_link().size() == 0)        // final round; nothing needs to be sent
       return;


    // step 2: enqueue
    for (int i = 0; i < rp.out_link().size(); ++i)    // redundant since size should equal to 1
    {
        // only send to root of group, but not self
        if (rp.out_link().target(i).gid != rp.gid())
        {
	  rp.enqueue(rp.out_link().target(i), b->state);
        }
        //else
            //fmt::print(stderr, "[{}:{}] Skipping sending to self\n", rp.gid(), round);
    }
}



// diy::decompose needs to have a function defined to create a block
// here, it is wrapped in an object to add blocks with an overloaded () operator
// it could have also been written as a standalone function

template <typename Bnds, typename Link, typename B>
struct ConfigBlockAdder
{
  typedef Link link_type;
  typedef Bnds bounds_type;
  typedef B block_type;

  ConfigBlockAdder(diy::Master& master_):
    master(master_) { }

  // this is the function that is needed for diy::decompose
  void  operator()(int gid,                // block global id
		   const bounds_type& core,     // block bounds without any ghost added
		   const bounds_type& bounds,   // block bounds including any ghost region added
		   const bounds_type& domain,   // global data bounds
		   const link_type& link)     // neighborhood
  {
       block_type* b = new block_type(core);
       link_type* l = new link_type(link);
       diy::Master& m = const_cast<diy::Master&>(master);
       m.add(gid, b, l); // add block to the master (mandatory)
  }

  diy::Master&  master;
};

#endif
