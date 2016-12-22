#ifndef _HiggsRecoil_hh_
#define _HiggsRecoil_hh_

#include <string>
#include <iostream>
#include <fstream>
#include <marlin/Processor.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TTree.h>
#include <TFile.h>

class TTree;

class HiggsRecoil  : public marlin::Processor
{
	public:

		Processor*  newProcessor() { return new HiggsRecoil ; }

		HiggsRecoil();

		~HiggsRecoil() {};

		void init();

		void processEvent( LCEvent * evtP );

		void end();

	protected:
		std::string _treeFileName;
		std::string _treeName;
		std::string _colName;
		std::string _colAdcVals;
		
		int _leptonID,_overwrite;
    int countmachine;

		unsigned int _eventNr;
		
		TFile *tree_file;
		TTree *eventlook;
		TTree *_outputTree;
		TTree *eutree;
		TTree *ctree;
		
		int num,pid;
		int ch;
		int recnum,_nMCP;

		float eppx,eppy,eppz,eped;
		float empx,empy,empz,emed;
		float uppx,uppy,uppz,uped;
		float umpx,umpy,umpz,umed;
		
		float recumpx,recumpy,recumpz,recumed;
		float recuppx,recuppy,recuppz,recuped;
		float jet1px,jet1py,jet1pz,jet1ed;
		float jet2px,jet2py,jet2pz,jet2ed;
		float y12,y23,y34,y45;
		
		int ppid;
		float ppx,ppy,ppz,ped;
				
		int ntal,nothch,nothnch,nchnum,chnum,npl,nce;
		int eventbin,eventfin;
		int fevent;
		int num_ep,num_em,num_up,num_um,num_othch;

		std::string _fileName;
		std::ostream *_output;
		std::string _histFileName;
};

#endif


             
