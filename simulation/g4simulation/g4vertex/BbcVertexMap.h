#ifndef __BBCVERTEXMAP_H__
#define __BBCVERTEXMAP_H__

#include "BbcVertex.h"

#include <phool/PHObject.h>
#include <map>
#include <iostream>

class BbcVertexMap : public PHObject {
  
public:

  typedef std::map<unsigned int, BbcVertex*>::const_iterator ConstIter;
  typedef std::map<unsigned int, BbcVertex*>::iterator            Iter;
  
  BbcVertexMap();
  virtual ~BbcVertexMap();

  void identify(std::ostream &os = std::cout) const;
  void Reset() {clear();}
  int  IsValid() const {return 1;}
  
  bool   empty()                   const {return _map.empty();}
  size_t size()                    const {return _map.size();}
  size_t count(unsigned int idkey) const {return _map.count(idkey);}
  void   clear();
  
  const BbcVertex* get(unsigned int idkey) const;
        BbcVertex* get(unsigned int idkey); 
        BbcVertex* insert(BbcVertex* vertex);
  size_t           erase(unsigned int idkey) {
    delete _map[idkey];
    return _map.erase(idkey);
  }

  ConstIter begin()                   const {return _map.begin();}
  ConstIter  find(unsigned int idkey) const {return _map.find(idkey);}
  ConstIter   end()                   const {return _map.end();}

  Iter begin()                   {return _map.begin();}
  Iter  find(unsigned int idkey) {return _map.find(idkey);}
  Iter   end()                   {return _map.end();}
  
private:
  std::map<unsigned int, BbcVertex*> _map;
    
  ClassDef(BbcVertexMap, 1);
};

#endif // __BBCVERTEXLIST_H__
