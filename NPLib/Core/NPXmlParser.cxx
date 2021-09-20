#include "NPXmlParser.h"
#include <stdio.h>
using namespace NPL;
using namespace NPL::XML;
////////////////////////////////////////////////////////////////////////////////
block::block(){};
block::~block(){};
int block::AsInt(std::string name){
  parameter p(name,"void");  
  auto it = m_parameters.find(p);
  if(it!=m_parameters.end()){
    std::string s = it->GetValue(); 
    return atoi(s.c_str());
    }
  else
    return -1000;
}

double block::AsDouble(std::string name){
  parameter p(name,"void");  
  auto it = m_parameters.find(p);
  if(it!=m_parameters.end()){
    std::string s = it->GetValue(); 
    return atof(s.c_str());
    }
  else
    return -1000;
}
std::string block::AsString(std::string name){
  parameter p(name,"void");  
  auto it = m_parameters.find(p);
  if(it!=m_parameters.end()){
    std::string s = it->GetValue(); 
    return s;
    }
  else
    return "void";
}
 

////////////////////////////////////////////////////////////////////////////////
parameter::parameter(){};
parameter::~parameter(){};

////////////////////////////////////////////////////////////////////////////////
Channel::Channel(){};
Channel::~Channel(){};

////////////////////////////////////////////////////////////////////////////////
void XmlParser::LoadFile(std::string file){
   // First create engine
   TXMLEngine* xml = new TXMLEngine;
   // Now try to parse xml file
   // Only file with restricted xml syntax are supported
   XMLDocPointer_t xmldoc = xml->ParseFile(file.c_str());
   if (xmldoc==0) {
      delete xml;
      return;
   }
   // take access to main node
   XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   // display recursively all nodes and subnodes
   LoadNode(xml, mainnode, 1);

   // Release memory before exit
   xml->FreeDoc(xmldoc);
   delete xml;
  }

//////////////////////////////////////////////////////////////////////////////////
void XmlParser::DisplayNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level){
   // this function display all accessible information about xml node and its children
   printf("%*c node: %s\n",level,' ', xml->GetNodeName(node));
   // display namespace
   XMLNsPointer_t ns = xml->GetNS(node);
   if (ns!=0)
      printf("%*c namespace: %s refer: %s\n",level+2,' ', xml->GetNSName(ns), xml->GetNSReference(ns));
   // display attributes
   XMLAttrPointer_t attr = xml->GetFirstAttr(node);
   while (attr!=0) {
       printf("%*c attr: %s value: %s\n",level+2,' ', xml->GetAttrName(attr), xml->GetAttrValue(attr));
       attr = xml->GetNextAttr(attr);
   }
   // display content (if exists)
   const char* content = xml->GetNodeContent(node);
   if (content!=0)
      printf("%*c cont: %s\n",level+2,' ', content);
   // display all child nodes
   XMLNodePointer_t child = xml->GetChild(node);
   while (child!=0) {
      DisplayNode(xml, child, level+2);
      child = xml->GetNext(child);
   }
}

//////////////////////////////////////////////////////////////////////////////////
void XmlParser::LoadNode(TXMLEngine* xml, XMLNodePointer_t node, Int_t level){
   // namespace
   XMLNsPointer_t ns = xml->GetNS(node);
   XMLNodePointer_t child = xml->GetChild(node);
   if(xml->GetNodeName(child)=="dataroot"){// main node
      std::cout <<" Loading XML file" << std::endl;
     }
   else{
    while(child!=0) {
      block b;
      // getting attribute:
      XMLNodePointer_t param = xml->GetChild(child);
      while (param!=0) {
        parameter p(xml->GetNodeName(param),xml->GetNodeContent(param));
        b.AddParameter(p); 
        param=xml->GetNext(param);
      }
      std::string name = xml->GetNodeName(child);
      b.SetName(name);
      m_blocks[name].push_back(b);
      child = xml->GetNext(child);
    }  
    std::cout << " -> XML file loaded for " <<m_blocks.size() << " detectors" << std::endl;
  }
}

//////////////////////////////////////////////////////////////////////////////////
std::vector<NPL::XML::block*> XmlParser::GetAllBlocksWithName(std::string name){
   std::vector<NPL::XML::block*> res;
   auto it=m_blocks.find(name);

    if(it!=m_blocks.end()){
      unsigned int size = it->second.size();
      for(unsigned int i = 0 ; i < size ; i++){
          res.push_back(&(it->second[i])); 
        }
      }
  
   return res;
  }
/////////////////////////////////////////////////////////////////////////////////
std::vector<std::string> XmlParser::GetAllBlocksName(){
  std::vector<std::string> res;
 for(auto it=m_blocks.begin(); it!= m_blocks.end();++it){
     res.push_back(it->first);
   } 
   return res;
  }
