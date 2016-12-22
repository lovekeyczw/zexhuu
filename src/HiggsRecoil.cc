#include <HiggsRecoil.hh>
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/LCFloatVec.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/MCParticleImpl.h>
#include <values.h>
#include <string>
#include <iostream>
#include <EVENT/LCFloatVec.h>
#include <EVENT/LCParameters.h>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TRandom.h>
#include <Rtypes.h>
#include <sstream>
#include <cmath>
#include <vector>
#include <TMath.h>
#include <fstream>
#include "TLorentzVector.h"
#include "IMPL/ReconstructedParticleImpl.h"
#include "UTIL/LCRelationNavigator.h"
#include "UTIL/PIDHandler.h"
#include "EVENT/LCFloatVec.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"

using namespace std;
using namespace fastjet;

const float sqrts = 250.0; 	//GeV
ofstream oft;
string fkk;
//ofstream ofp("total.dat");

int cont=0;

class pands
{
	public:
		int pos;
		ReconstructedParticle* p;
		
		pands(int post=-1,ReconstructedParticle* pp=NULL)
		{
			pos=post;
			p=pp;
		}
	};

/*void search(MCParticle *a)
{
	cont++;
	for(int i=0;i<cont;i++)
	{
		ofp<<"\t";
	}
	int nd= a->getDaughters().size();
	ofp<<a->getPDG()<<"\t";
	ofp<<a->getEnergy()<<"\t"<<a->getMomentum()[0]<<"\t"<<a->getMomentum()[1]<<"\t";
	ofp<<a->getParents().size()<<"\t";
	ofp<<nd<<endl;
	ofp<<a->isCreatedInSimulation()<<endl;
	for(int i=0;i<nd;i++)
	{
		ofp<<"\t";
		MCParticle *b=a->getDaughters()[i];
		search(b);
	}
	cont--;
}
*/

	int zhnum(LCCollection* col_MCP,int pdgpar)
	{
		int _nMCP = col_MCP->getNumberOfElements();
		int znum=-1;
		int hnum=-1;
		for (int i = 0; i < _nMCP; i++)
		{
			  MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(col_MCP->getElementAt(i));
				int tmpPID = a1_MCP->getPDG();
				int NDaughter = a1_MCP->getDaughters().size();
				MCParticle *daup, *daup2;
				
				
		/*		if(tmpPID==23&&NDaughter>=2)
					{
						//cout<<i<<endl;
					  int flagepc=0;
						int flagemc=0;
						MCParticle *ptmp=a1_MCP;
						
						
							for(int d=0;d<NDaughter;d++)
							{
								daup=ptmp->getDaughters()[d];
								int pdgtmp=daup->getPDG();
						//		cout<<pdgtmp<<endl;
							
								if(pdgtmp==pdgpar)
									{
										flagepc=1;
									}
								if(pdgtmp==-pdgpar)
									{
										flagemc=1;
									}
							}
						if(flagepc==1&&flagemc==1)
						{
							znum=i;
						}	
					}*/
						
					if(tmpPID==25&&NDaughter>=2)
					{
					//	cout<<"pid=25: "<<i<<endl;
					  int flagumc=0;
						int flagupc=0;
						MCParticle *ptmp=a1_MCP;
						
						
							for(int d=0;d<NDaughter;d++)
							{
								daup=ptmp->getDaughters()[d];
								int pdgtmp=daup->getPDG();
						//		cout<<pdgtmp<<endl;
							
								if(pdgtmp==pdgpar)
									{
										flagupc=1;
									}
								if(pdgtmp==-pdgpar)
									{
										flagumc=1;
									}
							}
						if(flagupc==1&&flagumc==1)
						{
							hnum=i;
						}	
					}
		}
					
					if(hnum!=-1)
						{
							return pdgpar;
						}
					else
						{
							return -1;
						}
	}
	
//ofstream ofm("test.dat");
//ofstream oftau("recma.dat");

	HiggsRecoil a_HiggsRecoil_instance;

	HiggsRecoil::HiggsRecoil()
		: Processor("HiggsRecoil"),
			_output(0)
			{
				_description = "Print MC Truth";

				_treeFileName = "MCTruth.root";
				registerProcessorParameter("TreeOutputFile",
				"The name of the file to which the ROOT tree will be written",
				_treeFileName,
				_treeFileName);

				_colName = "MCParticle";
				registerProcessorParameter("MCObjects",
				"The name of the PFOs",
				_colName,
				_colName);

				_treeName = "MCPart";
				registerProcessorParameter("TreeName",
				"The name of the ROOT tree",
				_treeName,
				_treeName);

			//	_leptonID = 13;k
				registerProcessorParameter("LeptonIDTag",
				"Lepton ID that will be used in this analysis.",
				_leptonID,
				_leptonID);

				_overwrite = 0;
				registerProcessorParameter("OverwriteFile",
				"If zero an already existing file will not be overwritten.",
				_overwrite,
				_overwrite);

			}


			void HiggsRecoil::init() {

				printParameters();
				
				tree_file = new TFile(_treeFileName.c_str(), "RECREATE");
				fkk=_treeFileName.replace(_treeFileName.find(".root"), 5, ".dat");

			/*		if (!tree_file->IsOpen()) {
						delete tree_file;
						tree_file = new TFile(_treeFileName.c_str(), "NEW");
					}*/

					_outputTree = new TTree(_treeName.c_str(), _treeName.c_str());
			    ctree=new TTree("ctree","ctree");
					eutree=new TTree("eutree","eutree");
					
					_outputTree->SetAutoSave(32 * 1024 * 1024);  // autosave every 32MB
					eutree->SetAutoSave(32 * 1024 * 1024);
				  ctree->SetAutoSave(32 * 1024 * 1024);
				  
				  // count how many event run
			  	ctree->Branch("nce", &nce, "nce/I");
				
					// store all particle
					_outputTree->Branch("ntal",  &ntal, "ntal/I");
					_outputTree->Branch("pid",   &ppid,  "ppid/I");
					_outputTree->Branch("ppx",   &ppx,  "ppx/F");
					_outputTree->Branch("ppy",   &ppy,  "ppy/F");
					_outputTree->Branch("ppz",   &ppz,  "ppz/F");
					_outputTree->Branch("ped",   &ped,  "ped/F");
					
					//some flag and num of photon and no photon
					eutree->Branch("recnum",    &recnum,    "recnum/I");
					eutree->Branch("fevent",    &fevent,    "fevent/I");
					eutree->Branch("npl",       &npl,       "npl/I"   );
					eutree->Branch("nothnch",   &nothnch,   "nothnch/I");
					eutree->Branch("numem",     &num_em,    "num_em/I");
					eutree->Branch("numep",     &num_ep,    "num_ep/I");
					eutree->Branch("numum",     &num_um,    "num_um/I");
					eutree->Branch("numup",     &num_up,    "num_up/I");
					eutree->Branch("numothch",  &num_othch, "num_othch/I");
					
					//store e+,e-,u+,u- 's Momentum and Energy
					eutree->Branch("jet1px",   &jet1px,  "jet1px/F");
					eutree->Branch("jet1py",   &jet1py,  "jet1py/F");
					eutree->Branch("jet1pz",   &jet1pz,  "jet1pz/F");
					eutree->Branch("jet1ed",   &jet1ed,  "jet1ed/F");
					eutree->Branch("jet2px",   &jet2px,  "jet2px/F");
					eutree->Branch("jet2py",   &jet2py,  "jet2py/F");
					eutree->Branch("jet2pz",   &jet2pz,  "jet2pz/F");
					eutree->Branch("jet2ed",   &jet2ed,  "jet2ed/F");
					eutree->Branch("recuppx",   &recuppx,  "recuppx/F");
					eutree->Branch("recuppy",   &recuppy,  "recuppy/F");
					eutree->Branch("recuppz",   &recuppz,  "recuppz/F");
					eutree->Branch("recuped",   &recuped,  "recuped/F");
					eutree->Branch("recumpx",   &recumpx,  "recumpx/F");
					eutree->Branch("recumpy",   &recumpy,  "recumpy/F");
					eutree->Branch("recumpz",   &recumpz,  "recumpz/F");
					eutree->Branch("recumed",   &recumed,  "recumed/F");
					eutree->Branch("y12",       &y12,      "y12/F");
					eutree->Branch("y23",       &y23,      "y23/F");
					eutree->Branch("y34",       &y34,      "y34/F");
					eutree->Branch("y45",       &y45,      "y45/F");
					//the begin and end of event
					eutree->Branch("eventbin",  &eventbin, "eventbin/I");
					eutree->Branch("eventfin",  &eventfin, "eventfin/I");
					
					eventbin=0;
					eventfin=0;
			
					countmachine = 0;
				//	cout<<"ini_ok"<<endl;
					ntal=0;
				}

				void HiggsRecoil::processEvent(LCEvent * evtP)
					{
						//cout<<"lll"<<endl;
						if (evtP)
							{
								try
								{
									ntal++;
								//	cout<<ntal<<endl;
									oft.open(fkk.c_str());  //the log of event num
									oft<<"event: "<<ntal<<endl;
									oft.close();
									
									LCCollection* col_MCP = evtP->getCollection("MCParticle");
									LCCollection* col_RecoP = evtP->getCollection( "ArborPFOs" );
								//	LCCollection* MrR = evtP->getCollection("RecoMCTruthLink");

									int _nMCP = col_MCP->getNumberOfElements();
									int _nRecoP = col_RecoP->getNumberOfElements();
									
							//		LCRelationNavigator *RCtoMC   = new LCRelationNavigator(MrR);

								//	cout<<"rec:  "<<_nRecoP<<endl;
								//	cout<<"mc:  "<<_nMCP<<endl;

									_eventNr = evtP->getEventNumber();
									//oft<<_eventNr<<endl;
									
						//			ofp<<"ntal:  "<<ntal<<endl;
									int tmpz=-1;

									for (int i = 0; i < _nMCP; i++)
									{
										  MCParticle *a1_MCP = dynamic_cast<EVENT::MCParticle *>(col_MCP->getElementAt(i));
											int tmpPID = a1_MCP->getPDG();
											if(tmpPID==25)
												{
													tmpz=abs(a1_MCP->getDaughters()[0]->getPDG());
												}
									}

									jet2px= 999;
									jet2py= 999;
									jet2pz= 999;
									jet2ed= 999;
									
									jet1px= 999;
									jet1py= 999;
									jet1pz= 999;
									jet1ed= 999;
									
									recuppx= 999;
									recuppy= 999;
									recuppz= 999;
									recuped= 999;
									
									recumpx= 999;
									recumpy= 999;
									recumpz= 999;
									recumed= 999;
									
									ppx=999;
									ppy=999;
									ppz=999;
									ped=999;
									pid=999;
								
									int zhn=-1;	

									
									nothch=0;
									nothnch=0;
									nchnum=0;
									chnum=0;
									npl=0;
									int nept=0;
									int nemt=0;
									int nupt=0;
									int numt=0;
									
									
									TLorentzVector ecms(0, 0, 0, 250);
									vector<pands> vum;
									vector<pands> vup;
									
									for (int i = 0; i < _nRecoP; i++)  //count the number of different particles
									{

										  ReconstructedParticle* rcp = dynamic_cast<ReconstructedParticle *>(col_RecoP->getElementAt(i));
										  int pid=rcp->getType();
										  int ch=rcp->getCharge();
										  
										  if(ch!=0)
										  	{
												  if(pid==13)
												  	{
												  		numt++;
												  		pands umtmp(i,rcp);
												  		vum.push_back(umtmp);
												  	}
												  else if(pid==-13)
												  	{
												  		nupt++;
												  		pands uptmp(i,rcp);
												  		vup.push_back(uptmp);
												  	}
												  else if(pid==11)
												  	{
												  		nemt++;
												  	}
												  else if(pid==-11)
												  	{
												  		nept++;
												  	}
												  else
												  	{
												  		nothch++;
												  	}
												}
											else
												{
													if(pid==22)
														{
															npl++;
														}
														else
															{
																nothnch++;
															}
												}
										  
									}
									
									/*
									for(int i=1;i<=16;i++)
									{
										zhn=zhnum(col_MCP,i);
										if(zhn!=-1)
											{
												tmpz=zhn;
											}
									}*/
									
									if(!(numt<=2&&nupt<=2&&nemt<=1&&nept<=1&&nothch==0)&&numt>0&&nupt>0)
										{
										//	cout<<"---------ok----------"<<endl;
										//	ofp<<ntal<<endl;
											recnum=ntal;
											
											float umpx,umpy,umpz,umed,uppx,uppy,uppz,uped;
											pands umt,upt;
											float min=999;
											float maxpt=-1;
											for(int umi=0;umi<vum.size();umi++)
											{
												umpx=vum[umi].p->getMomentum()[0];
												umpy=vum[umi].p->getMomentum()[1];
												umpz=vum[umi].p->getMomentum()[2];
												umed=vum[umi].p->getEnergy();
												
												TLorentzVector umtmp(umpx,umpy,umpz,umed);
												if(umtmp.Pt()>maxpt)
													{
														maxpt=umtmp.Pt();
														umt=vum[umi];
													}
											}
											
											maxpt=-1;
											
												for(int upi=0;upi<vup.size();upi++)
												{		
													uppx=vup[upi].p->getMomentum()[0];
													uppy=vup[upi].p->getMomentum()[1];
													uppz=vup[upi].p->getMomentum()[2];
													uped=vup[upi].p->getEnergy();
													
													TLorentzVector uptmp(uppx,uppy,uppz,uped);
													if(uptmp.Pt()>maxpt)
													{
														maxpt=uptmp.Pt();
														upt=vup[upi];
													}
												}
												
											recumpx=umt.p->getMomentum()[0];
											recumpy=umt.p->getMomentum()[1];
											recumpz=umt.p->getMomentum()[2];
											recumed=umt.p->getEnergy();
											
											recuppx=upt.p->getMomentum()[0];
											recuppy=upt.p->getMomentum()[1];
											recuppz=upt.p->getMomentum()[2];
											recuped=upt.p->getEnergy();
											
									/*		if(zhn>0)
												{
													fevent=recnum;
									//				cout<<"you: "<<ntal<<endl;
												}
											else
												{
													fevent=-1;
												}
										*/
										
								  		fevent=tmpz;
												
											vector<fastjet::PseudoJet> input_particles;	
												
											for (int i = 0; i < _nRecoP; i++)  
											{

												  ReconstructedParticle* rcp = dynamic_cast<ReconstructedParticle *>(col_RecoP->getElementAt(i));
												  int pid=rcp->getType();
												  int ch=rcp->getCharge();
												  
												  //count the number of different particles
												  
												  if(!(umt.pos==i||upt.pos==i))
												  	{
														  ppid=pid;
														  ppx=rcp->getMomentum()[0];
														  ppy=rcp->getMomentum()[1];
														  ppz=rcp->getMomentum()[2];
														  ped=rcp->getEnergy();
														  
														  fastjet::PseudoJet thisPtc(ppx,ppy,ppz,ped);
															input_particles.push_back(thisPtc);
														}
														
												  _outputTree->Fill();
												  eventfin++;
												  
												  
											 }
											
											if(input_particles.size()>=2)
												{
													JetDefinition jet_def(ee_kt_algorithm);
													ClusterSequence cs(input_particles, jet_def);
													vector<PseudoJet> jets = sorted_by_E(cs.exclusive_jets(2));
													
													y12=cs.exclusive_ymerge(1);
													y23=cs.exclusive_ymerge(2);
													y34=cs.exclusive_ymerge(3);
													y45=cs.exclusive_ymerge(4);
													
													jet1px=jets[0].px();
													jet1py=jets[0].py();
													jet1pz=jets[0].pz();
													jet1ed=jets[0].e();
													
													jet2px=jets[1].px();
													jet2py=jets[1].py();
													jet2pz=jets[1].pz();
													jet2ed=jets[1].e();
													
													num_em=nemt;
													num_ep=nept;
													num_um=numt;
													num_up=nupt;
													num_othch=nothch;

													 
													eutree->Fill();
													eventbin=eventfin+1;
												}
										}
										  
									
								}
									catch (lcio::DataNotAvailableException err) {}
							}
					}

							void HiggsRecoil::end()
								{
									//cout<<"total"<<endl;
									//ofp<<"tota//l::::::::::::::::::"<<endl;
								//	ofp<<totenv<<endl;
								//	ofp.close();
									oft.close();
								//	oftau.close();
								//	ofm.close();

									//cout << "countmachine=" << countmachine << endl;

									if (_outputTree) {
										nce=ntal;
										ctree->Fill();
									//	TFile *tree_file = _outputTree->GetCurrentFile(); //just in case we switched to a new file
										tree_file->Write();
//cout<<"tree"<<endl;
										delete tree_file;
									}
								}

							
