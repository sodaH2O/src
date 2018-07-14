// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 * This program was prepared by PSI.
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 *
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Particle/ParticleSpatialLayout.h"
#include "Particle/ParticleBConds.h"
#include "Index/NDIndex.h"
#include "Region/RegionLayout.h"
#include "Message/Communicate.h"
#include "Message/Message.h"
#include "Utility/IpplInfo.h"
#include "Utility/IpplStats.h"
#include "Utility/IpplMessageCounter.h"


// forward declarations
template <unsigned Dim> class FieldLayout;
class UserList;

/////////////////////////////////////////////////////////////////////
// constructor, from a FieldLayout
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::ParticleSpatialLayout(FieldLayout<Dim>& fl)
        : RLayout(fl)
{
    setup();			// perform necessary setup
}


/////////////////////////////////////////////////////////////////////
// constructor, from a FieldLayout
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::ParticleSpatialLayout(FieldLayout<Dim>& fl,
        Mesh& mesh)
        : RLayout(fl,mesh)
{
    setup();			// perform necessary setup
}


/////////////////////////////////////////////////////////////////////
// constructor, from a RegionLayout
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::ParticleSpatialLayout(const
        RegionLayout<T,Dim,Mesh>& rl) : RLayout(rl)
{
    setup();			// perform necessary setup
}


/////////////////////////////////////////////////////////////////////
// default constructor ... this does not initialize the RegionLayout,
// it will be instead initialized during the first update.
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::ParticleSpatialLayout() : RLayout()
{
    setup();			// perform necessary setup
}


/////////////////////////////////////////////////////////////////////
// perform common constructor tasks
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::setup()
{

    unsigned i;			// loop variable

	caching = false;

    // check ourselves in as a user of the RegionLayout
    RLayout.checkin(*this);

    // create storage for message pointers used in swapping particles
    unsigned N = Ippl::getNodes();
    SwapMsgList = new Message*[N];
    for (i = 0; i < Dim; ++i)
        SwapNodeList[i] = new bool[N];
    PutList = new std::vector<size_t>[N];

    // create storage for the number of particles on each node
    // and flag for empty node domain
    NodeCount = new size_t[N];
    EmptyNode = new bool[N];
    for (i = 0; i < N; ++i)
    {
        NodeCount[i] = 0;
        EmptyNode[i] = false;
    }

    // recalculate which nodes are our neighbors in each dimension
    if (RLayout.initialized())
        rebuild_neighbor_data();
}


/////////////////////////////////////////////////////////////////////
// destructor
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::~ParticleSpatialLayout()
{

    delete [] NodeCount;
    delete [] EmptyNode;
    delete [] SwapMsgList;
    for (unsigned int i=0; i < Dim; i++)
        delete [] (SwapNodeList[i]);
    delete [] PutList;

    // check ourselves out as a user of the RegionLayout
    RLayout.checkout(*this);
}




/////////////////////////////////////////////////////////////////////
// for each dimension, calculate where neighboring Vnodes and physical
// nodes are located, and create a list of this data.  This need only
// be updated when the layout changes.
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::rebuild_neighbor_data()
{
#ifdef OLD_SWAP_PARTICLES

    int d, j;			// loop variables
    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();
    Inform msg2all("rebuild_neighbor_data ", INFORM_ALL_NODES);

    // initialize the message list and initial node count
    for (d = 0; d < Dim; ++d)
    {
        NeighborNodes[d] = 0;
        for (j = 0; j < N; j++)
            SwapNodeList[d][j] = false;
    }

    // local Rnode iterators
    typename RegionLayout<T,Dim,Mesh>::iterator_iv localVN, endLocalVN = RLayout.end_iv();

    // check for no local Rnodes
    if (RLayout.begin_iv() == endLocalVN)
    {
        EmptyNode[myN] = true;
    }
    else
    {
        // we need to consider each spatial dimension separately
        for (d = 0; d < Dim; ++d)
        {
            // determine number of neighbor nodes in this spatial dimension
            for (localVN = RLayout.begin_iv(); localVN != endLocalVN; ++localVN)
            {
                // for each local Vnode, the domain to check equals the local Vnode dom
                // except for current dimension, which is the total Field domain size
                NDRegion<T,Dim> chkDom((*localVN).second->getDomain());
                chkDom[d] = (RLayout.getDomain())[d];
                // use the RegionLayout to find all remote Vnodes which touch the
                // domain being checked here
                typename RegionLayout<T,Dim,Mesh>::touch_range_dv touchingVN =
                    RLayout.touch_range_rdv(chkDom);
                typename RegionLayout<T,Dim,Mesh>::touch_iterator_dv tVN = touchingVN.first;
                for ( ; tVN != touchingVN.second; ++tVN)
                {
                    // note that we need to send a message to the node which contains
                    // this remote Vnode
                    unsigned vn = ((*tVN).second)->getNode();
                    if ( ! SwapNodeList[d][vn] )
                    {
                        SwapNodeList[d][vn] = true;
                        NeighborNodes[d]++;
                    }
                }
            }
        }
    }

    // there is extra work to do if there are multipple nodes, to distribute
    // the particle layout data to all nodes
    if (N > 1)
    {
        // At this point, we can send our domain status to node 0, and
        // receive back the complete domain status.
        int tag1 = Ippl::Comm->next_tag(P_SPATIAL_LAYOUT_TAG, P_LAYOUT_CYCLE);
        int tag2 = Ippl::Comm->next_tag(P_SPATIAL_RETURN_TAG, P_LAYOUT_CYCLE);
        Message *msg1, *msg2;
        if (myN != 0)
        {
            msg1 = new Message;
            // put local domain status in the message
            msg1->put(EmptyNode[myN]);
            // send this info to node 0
            int node = 0;
            Ippl::Comm->send(msg1, node, tag1);

            // receive back the complete domain status
            msg2 = Ippl::Comm->receive_block(node, tag2);
            msg2->get(EmptyNode);
            delete msg2;
        }
        else  			// do update tasks particular to node 0
        {
            // receive messages from other nodes describing their domain status
            int notrecvd = N - 1;	// do not need to receive from node 0
            while (notrecvd > 0)
            {
                // receive a message from another node.  After recv, node == sender.
                int node = Communicate::COMM_ANY_NODE;
                msg1 = Ippl::Comm->receive_block(node, tag1);
                msg1->get(EmptyNode[node]);
                delete msg1;
                notrecvd--;
            }

            // send info back to all the client nodes
            msg2 = new Message;
            msg2->put(EmptyNode, EmptyNode+N);
            Ippl::Comm->broadcast_others(msg2, tag2);
        }
    }

    // fix-up for empty nodes
    if (EmptyNode[myN])
    {
        // node with empty domain treats all non-empty nodes as neighbor
        // nodes along dimension 0
        for (j = 0; j < N; ++j)
        {
            if (!EmptyNode[j])
            {
                SwapNodeList[0][j] = true;
                NeighborNodes[0]++;
            }
        }
    }
    return;
#endif
}


/////////////////////////////////////////////////////////////////////
// Update the location and indices of all atoms in the given IpplParticleBase
// object.  This handles swapping particles among processors if
// needed, and handles create and destroy requests.  When complete,
// all nodes have correct layout information.
template <class T, unsigned Dim, class Mesh, class CachingPolicy>
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::update(
    IpplParticleBase< ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy> >& PData,
    const ParticleAttrib<char>* canSwap)
{

    unsigned N = Ippl::getNodes();
    unsigned myN = Ippl::myNode();
    size_t LocalNum   = PData.getLocalNum();
    size_t DestroyNum = PData.getDestroyNum();
    size_t TotalNum;
    int node;

    //Inform dbgmsg("SpatialLayout::update", INFORM_ALL_NODES);
    //dbgmsg << "At start of update:" << endl;
    //PData.printDebug(dbgmsg);

    // delete particles in destroy list, update local num
    PData.performDestroy();
    LocalNum -= DestroyNum;

    // set up our layout, if not already done ... we could also do this if
    // we needed to expand our spatial region.
    if ( ! RLayout.initialized())
        rebuild_layout(LocalNum,PData);

    // apply boundary conditions to the particle positions
    if (this->getUpdateFlag(ParticleLayout<T,Dim>::BCONDS))
        this->apply_bconds(LocalNum, PData.R, this->getBConds(), RLayout.getDomain());

    if(caching)
        this->updateCacheInformation(*this);

    // swap particles, if necessary
    if (N > 1 && this->getUpdateFlag(ParticleLayout<T,Dim>::SWAP))
    {
        // Now we can swap particles that have moved outside the region of
        // local field space.  This is done in several passes, one for each
        // spatial dimension.  The NodeCount values are updated by this routine.

        static IpplMessageCounterRegion swapCounter("Swap Particles");
        swapCounter.begin();
		//MPI_Pcontrol( 1,"swap_particles");

#ifndef OLD_SWAP_PARTICLES
        if (canSwap==0)
            LocalNum = new_swap_particles(LocalNum, PData);
        else
            LocalNum = new_swap_particles(LocalNum, PData, *canSwap);
#else
        if (canSwap==0)
            LocalNum = swap_particles(LocalNum, PData);
        else
            LocalNum = swap_particles(LocalNum, PData, *canSwap);
#endif

		//MPI_Pcontrol(-1,"swap_particles");
        swapCounter.end();
    }

    // Save how many local particles we have.
    TotalNum = NodeCount[myN] = LocalNum;

    // there is extra work to do if there are multipple nodes, to distribute
    // the particle layout data to all nodes
    if (N > 1)
    {
        // At this point, we can send our particle count updates to node 0, and
        // receive back the particle layout.
        int tag1 = Ippl::Comm->next_tag(P_SPATIAL_LAYOUT_TAG, P_LAYOUT_CYCLE);
        int tag2 = Ippl::Comm->next_tag(P_SPATIAL_RETURN_TAG, P_LAYOUT_CYCLE);
        if (myN != 0)
        {
            Message *msg = new Message;

            // put local particle count in the message
            msg->put(LocalNum);
            // send this info to node 0
            Ippl::Comm->send(msg, 0, tag1);

            // receive back the number of particles on each node
            node = 0;
            Message* recmsg = Ippl::Comm->receive_block(node, tag2);
            recmsg->get(NodeCount);
            recmsg->get(TotalNum);
            delete recmsg;
        }
        else  			// do update tasks particular to node 0
        {
            // receive messages from other nodes describing what they have
            int notrecvd = N - 1;	// do not need to receive from node 0
            TotalNum = LocalNum;
            while (notrecvd > 0)
            {
                // receive a message from another node.  After recv, node == sender.
                node = Communicate::COMM_ANY_NODE;
                Message *recmsg = Ippl::Comm->receive_block(node, tag1);
                size_t remNodeCount = 0;
                recmsg->get(remNodeCount);
                delete recmsg;
                notrecvd--;

                // update values based on data from remote node
                TotalNum += remNodeCount;
                NodeCount[node] = remNodeCount;
            }

            // send info back to all the client nodes
            Message *msg = new Message;
            msg->put(NodeCount, NodeCount + N);
            msg->put(TotalNum);
            Ippl::Comm->broadcast_others(msg, tag2);
        }
    }

    // update our particle number counts
    PData.setTotalNum(TotalNum);	// set the total atom count
    PData.setLocalNum(LocalNum);	// set the number of local atoms
	
	if(caching)
            this->updateGhostParticles(PData, *this);

    //dbgmsg << endl << "At end of update:" << endl;
    //PData.printDebug(dbgmsg);
}





/////////////////////////////////////////////////////////////////////
// print it out
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
std::ostream& operator<<(std::ostream& out, const ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>& L)
{

    out << "ParticleSpatialLayout, with particle distribution:\n    ";
    for (unsigned int i=0; i < (unsigned int) Ippl::getNodes(); ++i)
        out << "Number of particles " << L.getNodeCount(i) << "  ";
    out << "\nSpatialLayout decomposition = " << L.getLayout();
    return out;
}


/////////////////////////////////////////////////////////////////////
// Print out information for debugging purposes.
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::printDebug(Inform& o)
{

    o << "PSpatial: distrib = ";
    for (int i=0; i < Ippl::getNodes(); ++i)
        o << NodeCount[i] << "  ";
}


//////////////////////////////////////////////////////////////////////
// Repartition onto a new layout, if the layout changes ... this is a
// virtual function called by a UserList, as opposed to the RepartitionLayout
// function used by the particle load balancing mechanism.
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::Repartition(UserList* userlist)
{
    // perform actions to restructure our data due to a change in the
    // RegionLayout
    if (userlist->getUserListID() == RLayout.get_Id())
    {
        //Inform dbgmsg("ParticleSpatialLayout::Repartition", INFORM_ALL_NODES);
        //dbgmsg << "Repartitioning due to change in RegionLayout ..." << endl;

        // recalculate which nodes are our neighbors in each dimension
        rebuild_neighbor_data();
    }
}


//////////////////////////////////////////////////////////////////////
// Tell the subclass that the FieldLayoutBase is being deleted, so
// don't use it anymore
template < class T, unsigned Dim, class Mesh, class CachingPolicy >
void ParticleSpatialLayout<T,Dim,Mesh,CachingPolicy>::notifyUserOfDelete(UserList*)
{

    // really, nothing to do, since the RegionLayout we use only gets
    // deleted when we are deleted ourselves.
    return;
}


/***************************************************************************
 * $RCSfile: ParticleSpatialLayout.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:29 $
 * IPPL_VERSION_ID: $Id: ParticleSpatialLayout.cpp,v 1.1.1.1 2003/01/23 07:40:29 adelmann Exp $
 ***************************************************************************/
