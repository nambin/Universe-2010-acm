#include <stdafx.h>
#include <Illumination/IlluminationHeader/Primitive/KdTreeAccel.h>

using namespace illuminationspace;

void KdAccelNode::initLeaf(int *primNums, int np, MailboxPrim *mailboxPrims, MemoryArena &arena) {
	nPrims = np << 2;
	flags |= 3;

	if(np == 0)			onePrimitive = NULL;
	else if(np == 1)	onePrimitive = &mailboxPrims[primNums[0]];
	else {
		primitives = (MailboxPrim **)arena.Alloc(np * sizeof(MailboxPrim*));
		for(int i=0;i<np;++i)
			primitives[i] = &mailboxPrims[primNums[i]];
	}
}