#include "sys/sysinfo.h"
#include"sys/time.h"
#include<dirent.h>
#include <sys/mman.h> 
#include<sys/types.h>
#include<sys/stat.h>
#include<stdlib.h> 
#include<memory.h> 
#include<iostream> 
#include<string.h>
#include <fcntl.h>
#include<fstream>
#include <map>
#include <vector>
#include <pthread.h>
#include <unistd.h>
#include <error.h>
 
#define MAX_INT_DIG 32
#define CHR_NUM 24


using namespace  std;



struct InitBag
{
	string strFile;
	int *memeStart;
	InitBag()
	{
		strFile = "";
		memeStart = NULL;
	}
};


struct ChainBag
{
	bool bReverse;
	//bool bHasSun;
	int  nSrcStart;
	int  nSrcEnd;
	int  nDesStart;
	int  nDesEnd;
	int  nSpanChr;    //第二版中存储 跨染色体 的染色体
	string strDesChr; // stand chains 标记跨染色提标记 add by lyj 11.15
	//nSpanChr
	ChainBag()
	{
		bReverse = false;
		//bHasSun = false;
		nSrcStart = -1;
		nSrcEnd = -1;
		nDesStart =-1;
		nDesEnd = -1;
		nSpanChr = 0;
		strDesChr = "";
	}
	~ChainBag()
	{
		bReverse = false;
		//bHasSun = false;
		nSrcStart = -1;
		nSrcEnd = -1;
		nDesStart =-1;
		nDesEnd = -1;
		nSpanChr = 0;
		strDesChr = "";
	}
};

struct  CmySection
{
	std::map<int,int> map_Chain_1V1;
	std::map<int,ChainBag> map_Chain;
	CmySection()
	{
		map_Chain_1V1.clear();
		map_Chain.clear();
	}
	~CmySection()
	{
		map_Chain_1V1.clear();
		map_Chain.clear();
	}
};

typedef std::map<string,CmySection>::iterator mapChr_Chains_iterator;
	
typedef std::map<string,std::map<int,ChainBag> >::iterator ChrChain_iterator;
typedef std::map <int,int>::iterator Map_IntInt_iterator;
typedef std::map <int,ChainBag>::iterator Map_Chain_iterator;
typedef pair<map<int,ChainBag>::iterator,bool> MapRet;

//typedef std::map <string,SNPBag>::iterator Map_Rs38_iterator;
	

//typedef  pair<map<string,SNPBag>::iterator,bool> Map_RS_Ret;



/*****************************************************************************
** Function:     to_string(int _Val)
**  扩展string <int> // 原因是<string> 库里面没有
** Create Date:  2018.7.17
** Modify Time:  2018.7.17
** Author:        LYJ
** Version:      1.0
*******************************************************************************/

string to_string(int nVal)
{	
	char cBuf[2 * MAX_INT_DIG];
        snprintf(cBuf, sizeof(cBuf), "%d", nVal);
	return (string(cBuf));
}

int stoi(string str)
{
	//cout << str << endl;
	int num = atoi( str.c_str());
	return num;
}



/*****************************************************************************
** Function:    is_xxx_exist(string sFilename)
** 判断文件或者路径是否存在
** Create Date:  2018.7.17
** Modify Time:  2018.7.17
** Author:        LYJ
** Version:       1.0
*******************************************************************************/

 bool is_file_exist(string sFilename)
 {
 	if( 0 != access(sFilename.c_str(),F_OK))
 		return false;
  return true;
 }
 
 
bool is_dir_exist(string sDir)
 {
 	
 	  DIR *dirptr=opendir(sDir.c_str());
 	  if( NULL == dirptr)
 	  	return false;
 	  
 	  closedir(dirptr);
 	  return true;
 }
 
 


/*****************************************************************************
** Function:   GetMakeupChain(string chianFile,std::map<int,ChainBag> &map_Chain)
** chianFile ：参考对照关系文件
** Create Date:  2018.8.20
** Modify Time:  
** Author:        LYJ
** Version:      1.0
*******************************************************************************/

int GetMakeupChain_119(string &strResult,std::map<int,ChainBag> &map_Chain_D)
{
	std::map<int,ChainBag> map_Chain;
	int np1 = 0;
	int np2 = 0;
	ChainBag c_bag;
	int nPS = 0;
	int nPE = 0;
	int nCount  = 0 ;
	string str_Line = "";
  int nLen = strResult.length();
  //cout << nLen << endl;

	while(nPE != -1)  
	{  
		np1 = 0;
		np2 = 0;
		nPE = strResult.find('\n',nPS);
		if(nPE == -1 ) // 最后一条
		{
			if(c_bag.bReverse)
			{
		     int nt = -c_bag.nSrcStart;  c_bag.nSrcStart = -c_bag.nSrcEnd; c_bag.nSrcEnd = nt;			
	    }
	    MapRet ret = map_Chain.insert(std::make_pair(c_bag.nSrcStart,c_bag)); // 最后一条
	    if (!ret.second){ cout<<"insert error"<<endl;}
			continue;
		}
		str_Line = strResult.substr(nPS,nPE-nPS);
		nPS = nPE + 1;
		np2 = str_Line.find(" ",np1);
		if (np2 == -1)
			continue;
		nCount ++;
	
		string srcStart = str_Line.substr(np1,np2-np1);
		np1 = np2 +1;
		int nStart = stoi(srcStart);

		/////////////////////nSEND//////////////////////////
		np2 = str_Line.find(" ",np1);
		string srcEnd = str_Line.substr(np1,np2-np1);
		np1 = np2 +1;
		int nSEND = stoi(srcEnd);
		/////////////////////strand//////////////////////////
		np2 = str_Line.find(" ",np1);
		string sstrand = str_Line.substr(np1,np2-np1);
		np1 = np2 +1;
		
		if(sstrand == "-"){
			  nStart = -1 *nStart;
			  nSEND = -1*nSEND;
		}
		/////////////////////sDesStart//////////////////////////
		np2 = str_Line.find(" ",np1);
		string sDesStart = str_Line.substr(np1,np2-np1);
		np1 = np2 +1;
		int nES = stoi(sDesStart);
		//////////////////////sDesEnd/////////////////////////
		np2 = str_Line.find(" ",np1);
		string sDesEnd = str_Line.substr(np1,np2-np1);
		np1 = np2 +1;
		int nED = stoi(sDesEnd);
		
	  if ((nStart == c_bag.nSrcEnd +1)&&(nES == c_bag.nDesEnd +1) )// 说明是连起来的
		{
			c_bag.nSrcEnd = nSEND;
			c_bag.nDesEnd = nED;
		}
		else
		{
			if (c_bag.nSrcStart != -1)
			{
				if(c_bag.bReverse) // 交换
					{
						int nt = -c_bag.nSrcStart;
						c_bag.nSrcStart = -c_bag.nSrcEnd;
						c_bag.nSrcEnd = nt;
					}
				MapRet ret = map_Chain.insert(std::make_pair(c_bag.nSrcStart,c_bag));
				if (!ret.second){ cout<<"insert error"<<endl;}
			}
			
			c_bag.nSrcStart = nStart;
			c_bag.nSrcEnd = nSEND;
			if (sstrand == "-")
				c_bag.bReverse = true;
			else
				c_bag.bReverse = false;

			c_bag.nDesStart = nES;
			c_bag.nDesEnd = nED;
		}	
	}  //end while  挑出所有的对应关系 并做了排序
	
	//////////////////////处理反向连接将小区间连接起来//
	ChainBag c_bag2;
	bool bhas = false;
	std::map<int,ChainBag>::iterator it = map_Chain.begin();
	for(;it != map_Chain.end(); it++)
	{
		if(it->second.bReverse && !bhas) //第一包初始化
		{
			c_bag2.bReverse = true;
			c_bag2.nSrcStart = it->second.nSrcStart;
			c_bag2.nSrcEnd = it->second.nSrcEnd;
			c_bag2.nDesStart = it->second.nDesStart;
			c_bag2.nDesEnd = it->second.nDesEnd;
			bhas = true;
			continue;
		}
		else if(it->second.bReverse)
		{ 
			//(nStart == c_bag.nSrcEnd +1)&&(nES == c_bag.nDesEnd +1)
			if( (it->second.nSrcStart == c_bag2.nSrcEnd +1) && (it->second.nDesEnd +1 == c_bag2.nDesStart)) // 连起来的
				{
					c_bag2.nSrcEnd = it->second.nSrcEnd;
					c_bag2.nDesStart = it->second.nDesStart;
					continue;
				}
			else //不连续则插入
				{
					int nT = c_bag2.nDesStart;	c_bag2.nDesStart = c_bag2.nDesEnd;c_bag2.nDesEnd = nT; //交换
					map_Chain_D.insert(make_pair(c_bag2.nSrcStart,c_bag2));
					////////////////////////////////////////////////////
					c_bag2.bReverse = true;
			    c_bag2.nSrcStart = it->second.nSrcStart;
			    c_bag2.nSrcEnd = it->second.nSrcEnd;
			    c_bag2.nDesStart = it->second.nDesStart;
			    c_bag2.nDesEnd = it->second.nDesEnd;
					continue;
				}
		}
		else
			{
				map_Chain_D.insert(make_pair(it->second.nSrcStart,it->second));
			}
	
	}
	
	if(c_bag2.nDesStart != -1)
		{
				nPE = c_bag2.nDesStart;	c_bag2.nDesStart = c_bag2.nDesEnd;c_bag2.nDesEnd = nPE; //交换
	      map_Chain_D.insert(make_pair(c_bag2.nSrcStart,c_bag2));
		}


	return nCount;
}

/*****************************************************************************
** Function:      bWritResult
** 将CHR VCF文件 提取到的信息重新写入
** Create Date:  2018.7.04
** Modify Time:  2018.7.04
** Author:        LYJ
** Version:      1.0
*******************************************************************************/
bool bWritResult(FILE * &pFile,string &strInfo)
{
	int nLen = strInfo.length();
	int nWriten= 0;
	nWriten = fwrite (strInfo.c_str(), 1, strInfo.length(), pFile);

	if (nLen != nWriten)
	{
		printf("write error % d     实际写入长度 ： %d \n",nLen ,nWriten);
		return false;
	}
	return true;
}


int Repalce_char(string &str,const char &old_value,const string& new_value)
{
	int nRe = 0;
	for(int pos=0;pos!=-1;pos+=new_value.length())
	{
		if((pos=str.find(old_value,pos))!=-1)
		{
			//str.replace(pos,old_value.length(),new_value);
			str.replace(pos,1,new_value);
			nRe++;
		}
		else 
			break;
	}
	return nRe;
}

int Repalce_char(string &str,const string &old_value,const string& new_value)
{
	int nRe = 0;
	for(int pos=0;pos!=-1;pos+=new_value.length())
	{
		if((pos=str.find(old_value,pos))!=-1)
		{
			str.replace(pos,old_value.length(),new_value);
			nRe++;
		}
		else 
			break;
	}
	return nRe;
}

bool  getCurrentPath(string &sPath )
{
	char current_path[1024]; 
  int cnt = readlink("/proc/self/exe",current_path, 1024); 
  if (cnt < 0 || cnt >= 1024) 
  { 
    printf("***Error***\n"); 
     return false;
  } 
  int i; 
  for (i = cnt; i >=0; --i) 
  { 
    if (current_path[i] == '/') 
    { 
        current_path[i] = '\0'; 
        break; 
    }
  } 
  sPath = current_path ;
  return true;
}

  
int  GetChrSum(string sChr,int &nS19,int &nS38)
{
	 if(sChr == "chrX"){
	 		nS19 = 155270560;
      nS38 = 156040895;
      return 0;
	 	}
   if(sChr == "chrY"){
    		nS19 = 59373566;
        nS38 = 57227415;
        return 0;
    }
   if(sChr == "chr1"){
    		nS19 = 249250621;
        nS38 = 248956422;
        return 0;
    }
    if(sChr == "chr2"){
    		 nS19 = 243199373;
         nS38 = 242193529;
        return 0;
    }
    if(sChr == "chr3"){
    		 nS19 = 198022430;
         nS38 = 198295559;
        return 0;
    }
   if(sChr == "chr4"){
    		nS19 = 191154276;
        nS38 = 190214555;
        return 0;
    }
   if(sChr == "chr5"){
    		nS19 = 180915260;
        nS38 = 181538259;
        return 0;
   }
   if(sChr == "chr6"){
    		nS19 = 171115067;
        nS38 = 170805979;
        return 0;
   }
  if(sChr == "chr7"){
    		nS19 = 159138663;
        nS38 = 159345973;
        return 0;
   }
  if(sChr == "chr8"){
    		nS19 = 146364022;
        nS38 = 145138636;
        return 0;
   }
  if(sChr == "chr9"){
    		nS19 = 141213431;
        nS38 = 138394717;
        return 0;
   }
  if(sChr == "chr10"){
    		nS19 = 135534747;
        nS38 = 133797422;
        return 0;
   }
  if(sChr == "chr11"){
    		 nS19 = 135006516;
        nS38 = 135086622;
        return 0;
   }
  if(sChr == "chr12"){
    		nS19 = 133851895;
        nS38 = 133275309;
        return 0;
   }
  if(sChr == "chr13"){
    		nS19 = 115169878;
        nS38 = 114364328;
        return 0;
   }
  if(sChr == "chr14"){
    		nS19 = 107349540;
        nS38 = 107043718;
        return 0;
   }
  if(sChr == "chr15"){
    		nS19 = 102531392;
        nS38 = 101991189;
        return 0;
  }
  if(sChr == "chr16"){
    		nS19 = 90354753;
        nS38 = 90338345;
        return 0;
  }
  if(sChr == "chr17"){
    		nS19 = 81195210;
        nS38 = 83257441;
        return 0;
  }
  if(sChr == "chr18"){
    		nS19 = 78077248;
        nS38 = 80373285;
        return 0;
  }
  if(sChr == "chr19"){
    		 nS19 = 59128983;
         nS38 = 58617616;
        return 0;
  }
  if(sChr == "chr20"){
    		nS19 = 63025520;
        nS38 = 64444167;
        return 0;
  }
  if(sChr == "chr21"){
    		nS19 = 48129895;
        nS38 = 46709983;
        return 0;
  }
  if(sChr == "chr22"){
    		nS19 = 51304566;
        nS38 = 50818468;
        return 0;
   }
    
  return 1;
}


//通过对应关系组建标准链
int nchangechain2(std::map<int,ChainBag> &map_Chain,string &desf,string schr,int n19,int n38,int &nChain)
{
	//schr = "chr" + schr;
  int np1 = 0;
	int np2 = 0;
	int nCount  = 0 ;
	string str_Line = "";
	string sLastand = ""; //标记上一条的位置
	int nSrcS = 0;
	int nSrcE = 0;
	
	int ndesS = 0;  // h38
	int ndesE = 0;  // h38
	int nChainNum = 1;
	string sChainInfo = "";
	//string sSsum = to_string(n19);
	//string sDsum = to_string(n38);
	string sSrcS = ""; // h19 start
	string sSrcE = ""; // h19 End
	string sdesS = ""; // h38 start
	string sdesE = ""; // h38 End
	
	int nsrcLstE = 0;
	int ndeslstE = 0;
	
	string ssub = "";
	int nGap = 0;
	
  std::map<int,ChainBag>::iterator it = map_Chain.begin();
 for(;it != map_Chain.end();it++)
 {  
		int nStart = it->second.nSrcStart -1;
		int nSEND = it->second.nSrcEnd;
		string sstrand = "+";
		if(it->second.bReverse) {sstrand = "-";}
		int nES = it->second.nDesStart;
		int nED = it->second.nDesEnd;
		 nGap = nSEND - nStart;
		  
	  if(sLastand == "") 
		 {
        sSrcS = to_string(nStart);
        nsrcLstE = nSEND;
        ndeslstE =  nED;	
		  	
		  	if(sstrand == "+")
		  	{
		  		ndesS = nES -1;
		  		ndesE = nED;
		  	}
		  	else // 反向
		  	{ 
		  		ndesE = n38 - nES - 2;
		  		ndesS = ndesE + nGap;
		  	}
		    ssub = to_string(nGap);
		    sLastand = sstrand;
		    
		    if(sstrand == "-") // 直接输出
		    {
		    		 sChainInfo = "";
             sChainInfo = "chain 1 " + schr + " " + to_string(n19) + " + " + sSrcS + " " + to_string(nsrcLstE) + " " + schr ;
            if(ndesS > ndesE )
            {
             int nT = ndesS;
             ndesS = ndesE;
             ndesE = nT;
         }
		  	      //sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndeslstE) + " " + to_string(nChain++) +  "\n";
		  	      sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndesE) + " " + to_string(nChain++) +  "\n";
		  	      sChainInfo += ssub + "\n\n";
		  	     // cout << sChainInfo ;
		  	      desf += sChainInfo;
		  	     sLastand =  "";
		  	     continue;
		    }
		    
		    continue;			
		  }//第一条链
		  
		  if( sLastand == sstrand && sstrand == "+" ) // 同一条链
		  {
		  	int nG1 =  nStart - nsrcLstE ; // 偏移量
		  	int nG2 =  0;
        nG2 = nES - 1 - ndeslstE;
       
		  	if(nG2 < 0) // 换链
		  	 {
		  	 		  sChainInfo = "";
             	sChainInfo = "chain 1 " + schr + " " + to_string(n19) + " + " + sSrcS + " " + to_string(nsrcLstE) + " " + schr ;

		  	      //sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndeslstE) + " " + to_string(nChain++) +  "\n";
		  	      sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndesE) + " " + to_string(nChain++) +  "\n";
		  	      sChainInfo += ssub + "\n\n";
		  	      //cout << sChainInfo ;
		  	      desf += sChainInfo;
		  	      
		  	      sChainInfo = "";
		  	      nsrcLstE = 0;
		  	      ndeslstE = 0;
		  	      sLastand = "";
		  	      ssub = "";
		  	      {
		  	      	sLastand = sstrand;
                sSrcS = to_string(nStart);
                nsrcLstE = nSEND;
                ndeslstE =  nED;	
		  	
		  	       if(sstrand == "+")
		  		      {
		  			      ndesS = nES -1;
		  			      ndesE = nED;
		  		      }
		           ssub = to_string(nGap);
		            continue;	
		  	      }
		  	 }
		  	 ndesE = nED;
		  	 ndeslstE =  nED;
		  	 nsrcLstE = nSEND;
		     ssub += " " + to_string(nG1) + " " + to_string(nG2);
		     ssub += "\n" + to_string(nGap);	
		  }
		  else // 换链
		  {
		  	sChainInfo = "";
		  	sChainInfo = "chain 1 " + schr + " " + to_string(n19) + " + " + sSrcS + " " + to_string(nsrcLstE) + " " + schr ;
        if(ndesS > ndesE )
         {
             int nT = ndesS;
             ndesS = ndesE;
             ndesE = nT;
         }

		  	sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndesE) + " " + to_string(nChain++) +  "\n";
		  	sChainInfo += ssub + "\n\n";
		  	//cout << sChainInfo ;
		  	desf += sChainInfo;
		  	sChainInfo = "";
		  	nsrcLstE = 0;
		  	ndeslstE = 0;
		  	sLastand = "";
		  	
		  	{
		  	   sLastand = sstrand;
           sSrcS = to_string(nStart);
           nsrcLstE = nSEND;
           ndeslstE =  nED;	
		  	
		  	   if(sstrand == "+")
		  		  {
		  			   ndesS = nES -1;
		  			   ndesE = nED;
		  		  }
		  	    else // 反向
		  		  { 
		  			   ndesE = n38 - nES - 2;
		  			   ndesS = ndesE + nGap;
		  		  }
		        ssub = to_string(nGap);
		        continue;	
		   }
		  
		  }
		
 }  //end for
	
  if(ssub != "")
		{
			  sChainInfo = "chain 1 " + schr + " " + to_string(n19) + " + " + sSrcS + " " + to_string(nsrcLstE) + " " + schr ;
			  if(ndesS > ndesE )
         {
           int nT = ndesS;
           ndesS = ndesE;
           ndesE = nT;
          }
		  	sChainInfo += " " + to_string(n38) + " " + sLastand + " " + to_string(ndesS) + " " + to_string(ndesE) + " " + to_string(nChain++) +  "\n";
		  	sChainInfo += ssub + "\n\n";
		  	//cout << sChainInfo ;
		  	desf += sChainInfo;
		  	
		  	sChainInfo = "";
		}
	
	return 0;
	
}

static inline int nSplitStr2List(string &str,vector<string> &sL,const string splt)
{
	int npS = 0;
	int npE = 0;
	int nRe = 0;
	string strGet = "";
	sL.clear();
	for (;npE != -1;)
	{
		npE = str.find(splt,npS);
		if (npE == -1) // 最后的字串
		{
			strGet = str.substr(npS);//最后
			sL.push_back(strGet);
			nRe++ ;
			continue;
		}
		strGet = str.substr(npS,npE-npS);
		npS = npE +1;
		sL.push_back(strGet);
		nRe ++;
	}
	return nRe;
}



/*bool bRightFile(string strShed, string schr)
{
	bool bre = true;
	if (strShed == "00000000_X")
	{
		if (schr != "chr1")
			bre = false;
	}
	else if (strShed == "00000001_X")
	{
		if (schr != "chr10")
			bre = false;
	}
	else if (strShed == "00000002_X")
	{
		if (schr != "chr11")
			bre = false;
	}
	else if (strShed == "00000003_X")
	{
		if (schr != "chr12")
			bre = false;
	}
	else if (strShed == "00000004_X")
	{
		if (schr != "chr13")
			bre = false;
	}
	else if (strShed == "00000005_X")
	{
		if (schr != "chr14")
			bre = false;
	}
	else if (strShed == "00000006_X")
	{
		if (schr != "chr15")
			bre = false;
	}
	else if (strShed == "00000007_X")
	{
		if (schr != "chr16")
			bre = false;
	}
	else if (strShed == "00000008_X")
	{
		if (schr != "chr17")
			bre = false;
	}
	else if (strShed == "00000009_X")
	{
		if (schr != "chr18")
			bre = false;
	}
	else if (strShed == "00000010_X")
	{
		if (schr != "chr19")
			bre = false;
	}
	else if (strShed == "00000011_X")
	{
		if (schr != "chr2")
			bre = false;
	}
	else if (strShed == "00000012_X")
	{
		if (schr != "chr20")
			bre = false;
	}
	else if (strShed == "00000013_X")
	{
		if (schr != "chr21")
			bre = false;
	}
	else if (strShed == "00000014_X")
	{
		if (schr != "chr22")
			bre = false;
	}
	else if (strShed == "00000015_X")
	{
		if (schr != "chr3")
			bre = false;
	}
	else if (strShed == "00000016_X")
	{
		if (schr != "chr4")
			bre = false;
	}
	else if (strShed == "00000017_X")
	{
		if (schr != "chr5")
			bre = false;
	}
	else if (strShed == "00000018_X")
	{
		if (schr != "chr6")
			bre = false;
	}
	else if (strShed == "00000019_X")
	{
		if (schr != "chr7")
			bre = false;
	}
	else if (strShed == "00000020_X")
	{
		if (schr != "chr8")
			bre = false;
	}
	else if (strShed == "00000021_X")
	{
		if (schr != "chr9")
			bre = false;
	}
	else if (strShed == "00000022_X")
	{
		//if (schr != "MZ")
			bre = false;
	}
	else if (strShed == "00000023_X")
	{
		if (schr != "chrX")
			bre = false;
	}
	else if (strShed == "00000024_X")
	{
		if (schr != "chrY")
			bre = false;
	}
	else
	{
		//cout<<"find special :  "<<strShed << endl;
		bre = false;
	}

	return bre;

}*/

bool bRightFile2(string strShed, string schr)
{
	//cout << "strShed: " << strShed ;
	//cout << "chr: " << schr << endl;
	
	Repalce_char(strShed,"_X","");
	Repalce_char(schr,"chr","");
	int n1 = strShed.find_first_not_of("0");
	
	if(n1 < 0 )
		n1 = 1;
	else
		{
			strShed = strShed.substr(n1);
	    n1 = stoi(strShed)+1;
		}
  
  int n2 = -1;
  if(schr == "X")
  	n2 = 23;
  else if(schr == "Y")
  	n2 = 24;
  else
  	n2 = stoi(schr);
  	
  if(n1 == n2)
  	return true;
  
  return false;

}

int nCountChar(string &str,const char &c)
{
	int nRe = 0;
	for(int pos=0;pos!=-1;pos++)
	{
		if((pos=str.find(c,pos))!=-1)
			nRe++;
		else 
			break;
	}
	return nRe;
}


string   PickPairsF(string &sQUERY,string &sINDEX,string &sREFERENCE,string &sSTRAND,int nPosbase)//成功抽取前向Q插入
{
	string stRE = "";
	int nSs = 0;// source start
	int nSe = 0;// source end
	int nDs = 0;// Des  end
	int nDe = 0;// Des end



	int ns1= sQUERY.find_first_not_of(" ");
	int ns2 = sQUERY.find(" ",ns1);
	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	string strS = sQUERY.substr(ns1,ns2 -ns1);
	nSs = stoi(strS )+ nPosbase;

	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	
	string sSrcBase = sQUERY.substr(ns1,ns2 -ns1);


	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	if (ns2 == -1)// 特殊格式    QUERY:       5001  5000       REFERENCE:    2768454  2768453   
		return stRE;

	string strSe = sQUERY.substr(ns1,ns2-ns1);
	nSe = stoi(strSe)+nPosbase;

	ns1= sREFERENCE.find_first_not_of(" ");
	ns2 = sREFERENCE.find(" ",ns1);
	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string strDS = sREFERENCE.substr(ns1,ns2 -ns1);
	nDs = stoi(strDS);

	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string sDesBase = sREFERENCE.substr(ns1,ns2 -ns1);

	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string strdSe = sREFERENCE.substr(ns1,ns2-ns1);
	nDe = stoi(strdSe);


	int nBaselen  = sDesBase.length();
	if (sDesBase.length() !=sSrcBase.length())
	{
		cout <<"base error"<<endl;
		stRE = "";
		return stRE;
		//system("pause");
	}

	// you 全是"-"的情况

	if (sDesBase.length() == nCountChar(sDesBase,'-') || sSrcBase.length() == nCountChar(sSrcBase,'-')  )
		return stRE;

	if (sINDEX.find("-") == -1)//说明此段是对齐的
	{
		stRE = to_string(nSs) + " ";
		stRE += to_string(nSe) + " ";
		stRE += sSTRAND  + " ";
		stRE += to_string(nDs) + " ";
		stRE += to_string(nDe)+"\n";
	}
	else //说明此段不是对齐的
	{
		int nlen1 = sSrcBase.length();
		char s,d;
		int nSs1 = 0; //保存37号偏移
		int nDs1 = 0; //保存38号偏移
		int nSrcEND = 0;
		int nDesEND = 0;
		char cLS;
		char cLD;

		for (int i = 0;i<nlen1;i++)
		{
			s = sSrcBase[i];// debug
			d = sDesBase[i];

			if (s == '-'&& d != '-')// 37号上偏
			{
				nSs1--;
				nSrcEND = nSs + nSs1;
				nDesEND = nDs + nDs1;

				if (nDs !=nDesEND && nSrcEND == nSs && cLS != '-') //单对单
				{
					stRE += to_string(nSs) + " ";
					stRE += to_string(nSrcEND) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND-1) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else if (nDs ==nDesEND || nSrcEND == nSs) //连续的“-”
				{  
					
					if (nSrcEND == nSs-1)//说明从开头就是“-”
					{
            nSs = nSrcEND +1;
					   nDesEND += 1;
					}
					else				
					   nSs = nSrcEND;
					nDs = nDesEND;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else// 中间偏移
				{
					stRE += to_string(nSs) + " ";
					stRE += to_string(nSrcEND) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND-1) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}

			}
			else if (s!= '-'&& d == '-')//38号偏移
			{
				nDs1--;
				nSrcEND = nSs + nSs1;
				nDesEND = nDs + nDs1;

				if (nDs ==nDesEND && nSrcEND != nSs && cLD != '-') //单对单
				{

					stRE += to_string(nSs) + " ";
					stRE += to_string(nSrcEND-1) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else if (nDs ==nDesEND || nSrcEND == nSs) // 38号开头偏移
				{
					if (nDesEND == nDs-1)//说明从开头就是“-”
					{
						nDs = nDesEND +1;
						nSrcEND += 1;
					}
					else
						nDs = nDesEND;
					nSs = nSrcEND;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else// 中间偏移
				{
					stRE += to_string(nSs) + " ";
					stRE += to_string(nSrcEND-1) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}

			}

			nSs1++;
			nDs1++;
			cLS  = s;
			cLD = d;
		}//end for
		//最后一段
		nSrcEND = nSs + nSs1;
		nDesEND = nDs + nDs1;
		if (nSrcEND != nSs && nDesEND != nDs)
		{
			stRE += to_string(nSs) + " ";
			stRE += to_string(nSrcEND-1) + " ";
			stRE += sSTRAND  + " ";
			stRE += to_string(nDs) + " ";
			stRE += to_string(nDesEND-1) + "\n";
		}
		////////////////////验证关系是否正确///////////////////////
		if (nSrcEND - 1 != nSe || nDesEND-1 != nDe)
		{
			cout<<" Tail deals error :  nSrcEnd = " << nSrcEND -1 <<  "nSe = " << nSe <<endl;
			cout<<" Tail deals error :  nSrcEnd = " << nDesEND-1<<  "nSe = " << nDe <<endl;
			cout << stRE << endl;
			//system("pause");
			
		}
		///////////////////////////////////////////////////////////

	}

	return stRE;

}

string   PickPairsC(string &sQUERY,string &sINDEX,string &sREFERENCE,string &sSTRAND,int nPosbase)//成功抽取反向互补插入
{
	string stRE = "";
	int nSs = 0;// source start
	int nSe = 0;// source end
	int nDs = 0;// Des  end
	int nDe = 0;// Des end

	

	int ns1= sQUERY.find_first_not_of(" ");
	int ns2 = sQUERY.find(" ",ns1);
	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	string strS = sQUERY.substr(ns1,ns2 -ns1);
	nSs = stoi(strS )+ nPosbase;

	int Y = nSs;
	nSs = Y- nSs;


	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	string sSrcBase = sQUERY.substr(ns1,ns2 -ns1);


	ns1 = sQUERY.find_first_not_of(" ",ns2);
	ns2 = sQUERY.find(" ",ns1);
	if (ns2 == -1)// 特殊格式    QUERY:       5001  5000       REFERENCE:    2768454  2768453   
		return stRE;
		
	string strSe = sQUERY.substr(ns1,ns2-ns1);
	nSe = stoi(strSe)+nPosbase;

	nSe = Y- nSe;

	ns1= sREFERENCE.find_first_not_of(" ");
	ns2 = sREFERENCE.find(" ",ns1);
	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string strDS = sREFERENCE.substr(ns1,ns2 -ns1);
	nDs = stoi(strDS);

	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string sDesBase = sREFERENCE.substr(ns1,ns2 -ns1);

	ns1 = sREFERENCE.find_first_not_of(" ",ns2);
	ns2 = sREFERENCE.find(" ",ns1);
	string strdSe = sREFERENCE.substr(ns1,ns2-ns1);
	nDe = stoi(strdSe);

	int nBaselen  = sDesBase.length();

	if (sDesBase.length() !=sSrcBase.length())
	{
		cout <<"base error"<<endl;
		//system("pause");
	}

	// you 全是"-"的情况
	if (sDesBase.length() == nCountChar(sDesBase,'-') || sSrcBase.length() == nCountChar(sSrcBase,'-')  )
		return stRE;


	if (sINDEX.find("-") == -1)//说明此段是对齐的
	{
		stRE = to_string(Y-nSs) + " ";
		stRE += to_string(Y-nSe) + " ";
		stRE += sSTRAND  + " ";
		stRE += to_string(nDs) + " ";
		stRE += to_string(nDe)+ "\n";
	}
	else //说明此段不是对齐的
	{
		int nlen1 = sSrcBase.length();
		char s,d;
		int nSs1 = 0; //保存37号偏移
		int nDs1 = 0; //保存38号偏移
		int nSrcEND = 0;
		int nDesEND = 0;
		char cLS;
		char cLD;

		for (int i = 0;i<nlen1;i++)
		{
			s = sSrcBase[i];// debug
			d = sDesBase[i];

			if (s == '-'&& d != '-')// 37号上偏
			{
				nSs1--;
				nSrcEND = nSs + nSs1;
				nDesEND = nDs + nDs1;

				if (nDs !=nDesEND && nSrcEND == nSs && cLS != '-') //单对单
				{
					stRE += to_string(Y-nSs) + " ";
					stRE += to_string(Y-nSrcEND) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND-1) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else if (nDs ==nDesEND || nSrcEND == nSs) //连续的“-”
				{  

					if (nSrcEND == nSs-1)//说明从开头就是“-”
					{
						nSs = nSrcEND +1;
						nDesEND += 1;
					}
					else				
						nSs = nSrcEND;
					nDs = nDesEND;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else// 中间偏移
				{
					stRE += to_string(Y-nSs) + " ";
					stRE += to_string(Y-nSrcEND) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND-1) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}

			}
			else if (s!= '-'&& d == '-')//38号偏移
			{
				nDs1--;
				nSrcEND = nSs + nSs1;
				nDesEND = nDs + nDs1;

				if (nDs ==nDesEND && nSrcEND != nSs && cLD != '-') //单对单
				{

					stRE += to_string(Y-nSs) + " ";
					stRE += to_string(Y-(nSrcEND-1)) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else if (nDs ==nDesEND || nSrcEND == nSs) // 38号开头偏移
				{
					if (nDesEND == nDs-1)//说明从开头就是“-”
					{
						nDs = nDesEND +1;
						nSrcEND += 1;
					}
					else
						nDs = nDesEND;
					nSs = nSrcEND;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}
				else// 中间偏移
				{
					stRE += to_string(Y-nSs) + " ";
					stRE += to_string(Y-(nSrcEND-1)) + " ";
					stRE += sSTRAND  + " ";
					stRE += to_string(nDs) + " ";
					stRE += to_string(nDesEND) + "\n";
					nSs = nSrcEND+1;
					nDs = nDesEND+1;
					nSs1 = 0;
					nDs1 = 0;
					continue;
				}

			}

			nSs1++;
			nDs1++;
			cLS  = s;
			cLD = d;
		}//end for
		//最后一段
		nSrcEND = nSs + nSs1;
		nDesEND = nDs + nDs1;
		if (nSrcEND != nSs && nDesEND != nDs)
		{
			stRE += to_string(Y-nSs) + " ";
			stRE += to_string(Y-(nSrcEND-1)) + " ";
			stRE += sSTRAND  + " ";
			stRE += to_string(nDs) + " ";
			stRE += to_string(nDesEND-1) + "\n";
		}
		////////////////////验证关系是否正确///////////////////////
		if (   nSe != (nSrcEND-1) ||   nDe != (nDesEND-1))
		{
			cout<<" Tail deals error :  nSrcEnd = " << (nSrcEND-1) <<  "   nSe = " << nSe <<endl;
			cout<<" Tail deals error :  nDESEnd = " << nDesEND-1<<  "  nDe = " << nDe <<endl;
			cout << stRE << endl;
			//system("pause");

		}
		///////////////////////////////////////////////////////////

	}

	return stRE;
}


int filterShredout(string sfname,FILE * &pFile)
{
	string  str_Line = "";
	string str_Chain = "";
	static int nChain = 2;

	string  sfilePath;
	ifstream ifs;
	int nCout =0;
	sfilePath = sfname.substr(0,sfname.rfind('/')+1);// 存储路径

  
	ifs.open(sfname.c_str() );
	if (!ifs)
	{
		cout<<"Read File Error, please check file name is right: \n"<<endl;
		return -1;
	}
  cout << "start convert \n";
	bool bwrite = false;
	bool bStartline = false;
	string schr = "";// 标识当前的chr
	string slastchr = "";// 标识当前的chr
	string sQuality = "";
	string sStrand = "";//反向互补标识
	string sResult = "";
	int nSrcHead = 0;
	int nSrcEnd = 0;
	int nPosbase = 0;
	int nDesHead = 0;
	int nDesEnd = 0;
	//int nNum =0;
	int nS19 = 0;
  int nS38 = 0;
	while(!ifs.eof())  
	{  
		str_Line = "";
		getline(ifs,str_Line);
		if (str_Line.find("alignment") != -1 )
		{
			bwrite = false;
			bStartline = false;
			int npos = str_Line.find("S:");
			if (npos == -1)
				continue;
			npos +=2;
			sQuality = str_Line.substr(npos,str_Line.find(" ",npos)-npos);
			sQuality = "#" + sQuality +  "# "; // 加入# 区分
			

			string sChFile = "";
			npos = str_Line.find("_X",0);
			if (npos == -1)//  alignment 没有_x;
				continue;
			sChFile = str_Line.substr(0,npos+2);
			npos = sChFile.rfind(' ');
			if (npos == -1)
			{
				cout <<"pos == -1:  " <<str_Line << endl;
				continue;
			}
			sChFile = sChFile.substr(npos+1,sChFile.length()-npos);

			npos = str_Line.find("chr",0);
			if (npos == -1)
			{
				cout <<"find chr ==-1  :" <<str_Line << endl;
				continue;
			}
				
			string strch = str_Line.substr(npos+3,2);
			Repalce_char(strch," ","");
			Repalce_char(strch,"	","");

			if(!bRightFile2(sChFile,"chr"+strch)){continue;}
			
			//alignment:S:50 5000  00000000_X000010000 chr1        1     5000     10001     15000   F    5000 100.00 5000 248956422
			int nps = str_Line.find("_X")+1;
			if (nps != -1)
			{
				bStartline = true;//对齐行的汇总
				int np3 = str_Line.find(" ",nps);//spc
				int np4 = str_Line.find_first_not_of('0',nps+1);
				np3 = str_Line.find(' ',np4);
				//alignment:S:42 5000  00000008_X000000000 chr17        1     5000    150208    155207   F    5000 100.00 5000 83257441
				string strPosStart = str_Line.substr(np4,np3-np4);
				if (strPosStart == "")
					nPosbase = 0;
				else
					nPosbase = stoi(strPosStart);

				int np2 = str_Line.find(' ',npos);//chr 后的空格
				schr = str_Line.substr(npos,np2-npos);
				
				if(slastchr == "")
					slastchr = schr;
			  else if(slastchr != schr)
			  	{			  		
			  		std::map<int,ChainBag> map_Chain;
			  		GetMakeupChain_119(sResult,map_Chain);
			  		sResult = "";
            if( 0 != GetChrSum(slastchr,nS19,nS38))
            {
            		cout << "Get Error Chr: " << slastchr << endl;
            		return -1;
            }
            
            nchangechain2(map_Chain,sResult,slastchr,nS19,nS38,nChain);
            str_Chain += sResult;
            
            sResult = "";
			  		slastchr = schr;
			  		/////////////////////////////////////////
			  	}
				

				np3 = str_Line.find_first_not_of(' ',np2);
				np4 = str_Line.find(' ',np3+1);
				string sStart = str_Line.substr(np3,np4-np3);
				nSrcHead = stoi(sStart);

				np3 = str_Line.find_first_not_of(' ',np4);
				np4 = str_Line.find(' ',np3+1);
				string sEND = str_Line.substr(np3,np4-np3);
				nSrcEnd = stoi(sEND);

				np3 = str_Line.find_first_not_of(' ',np4);
				np4 = str_Line.find(' ',np3+1);
				string sDesStart = str_Line.substr(np3,np4-np3);
				nDesHead = stoi(sDesStart);

				np3 = str_Line.find_first_not_of(' ',np4);
				np4 = str_Line.find(' ',np3+1);
				string sDesEnd = str_Line.substr(np3,np4-np3);
				nDesEnd = stoi(sDesEnd);

				bwrite = true;
				np3 = str_Line.find_first_not_of(' ',np4);
				np4 = str_Line.find(' ',np3+1);
				sStrand = str_Line.substr(np3,np4-np3);
			}

			nCout ++;
			str_Line = "";
		}// end of alignment 行处理

		if ((bwrite&&str_Line.length()>3) )
		{
			string sQUERY = str_Line;
			string sINDEX = "";
			string sREFERENCE = "";

			getline(ifs,sINDEX);
			getline(ifs,sREFERENCE);

			if (sREFERENCE.find("REFERENCE") == -1){
				cout << "Format REFERENCE error: "<<sREFERENCE<<endl;
			}
			if (sStrand == "F"){
					string ss = "+";
					str_Line = PickPairsF(sQUERY,sINDEX,sREFERENCE,ss,nPosbase);// 通过3行对应的数据进行处理
				}
			else{
					string ss = "-";
					str_Line = PickPairsC(sQUERY,sINDEX,sREFERENCE,ss,nPosbase);// 通过3行对应的数据进行处理
				}
				
			if (str_Line != "")
			{
				sResult += str_Line;
				//cout << schr << "\t" << str_Line;
			}

		}

	}  //end while
	
	if(sResult != "")
	{
		 std::map<int,ChainBag> map_Chain;
		 GetMakeupChain_119(sResult,map_Chain);
		 sResult = ""; //清空结果
		//cout << "*"<<slastchr << "*\n";
		//cout << map_Chain.size()<< endl;
			  		
     if( 0 != GetChrSum(slastchr,nS19,nS38)){
      cout << "Get Error Chr: " << slastchr << endl;
       return -1;
     }
            
     nchangechain2(map_Chain,sResult,slastchr,nS19,nS38,nChain);
     str_Chain += sResult;
     sResult = "";
	}
	
	
	 if( false == bWritResult(pFile,str_Chain))
	 	return -1;
	

	return 0;
}

int main(int arg,char *args[])
{
	
	string srcfile = "";
	string strDesFile = "";
  FILE *pFile = NULL;
   if(arg <2)
	 	{
	 		 return 0;
	 	}

	///////////////////////////获取输入参数////////////////////////////////////////////////////
	for(int i=0; i<arg; i++)
	{
   	string sGet = args[i];
		if(i == 1)
		{  
			if(true ==  is_file_exist(sGet))
			 {
				srcfile = sGet;
			 }
			 else
			 	{
			 		cout << "Source file not exist " << sGet << endl;
			 		return -1;
			 	}
				
		}
	 if(i == 2)
		{
			if(true != is_file_exist(sGet))
				{
					strDesFile = sGet;
					if ((pFile = fopen(strDesFile.c_str(),"a+"))==NULL)
	          {
		           perror("Error occurs when creat file");
		           return false;
	          }
				} 
			else
			{
			 	cout << "Output file exist " << sGet << endl;
			 	return -1;
			}
		}

	}
 
   if(  -1 != filterShredout(srcfile,pFile)) 
   	{
   		if (fclose (pFile) != 0)
	      {
		      perror("Error occurs when close file"); 
		      return false;
	      }
	      cout << "work finished \n";
   	}
   	else
   		{
   			fclose (pFile);
   			remove(strDesFile.c_str());
   		}
    
  
  	return 0;
	 
}


