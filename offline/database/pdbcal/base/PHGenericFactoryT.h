#ifndef PDBCAL_BASE_PHGENERICFACTORYT_H
#define PDBCAL_BASE_PHGENERICFACTORYT_H

#include <iostream>
#include <map>
#include <string>

/** A Factory to build objects.
 *
 * This is more or less a (somewhat) de-templatized version of what 
 * can be found in Loki class library 
 * ( "Modern C++ Design" by Andrei Alexandrescu 
 * (Addison Wesley ISBN 0-201-70431-5) chapter 8 on Object Factories )
 *
 * The (partial) de-templatization is only due to the fact 
 * that there's problably no chance to get the original Loki .h file 
 * parsed by CINT in order to get a decent ROOT dictionary 
 * (that's a pity because it's a nice piece of code ;-) )
 *
 */

template
< 
  class AbstractProduct
>
class DefaultFactoryError
{
protected:
  AbstractProduct* OnUnknownType(const std::string& /*id*/) { return 0; }
};

template
<
  class AbstractProduct,
  template<class> class FactoryErrorPolicy = DefaultFactoryError
>
class PHGenericFactoryT : public FactoryErrorPolicy<AbstractProduct>
{
public:

  typedef AbstractProduct* (*ProductCreator)();
  typedef std::string IdentifierType;

  class ProductCreatorPair
  {
  public:
    ProductCreatorPair()
      : creator_(), name_() {}
    ProductCreatorPair(ProductCreator creator, const char* productname)
      : creator_(creator),name_(productname) {}

    ProductCreator creator() const { return creator_; }
    std::string productname() const { return name_; }

  private:
    ProductCreator creator_;
    std::string name_;
  };

public:

  /// The factory is a singleton.
  static PHGenericFactoryT<AbstractProduct,FactoryErrorPolicy>& instance()
  { 
    static PHGenericFactoryT<AbstractProduct,FactoryErrorPolicy> factory;
    return factory; 
  }
  
  /// Create an object identified by the string id.
  AbstractProduct* create(const char* id)
  {
    // the id parameter is a const char* instead of std::string
    // just because we want to be able to access this one from
    // ROOT prompt, and Cint seems to have some problems with
    // std::string.

    typename CreatorMap::const_iterator it = 
      fCreators.find(IdentifierType(id));

    if ( it != fCreators.end() )
      {
	// invoke the creation function
	return (it->second.creator())();
      }

    // Upon failure we delegate to the FactoryError class
    // that might simply returns 0 or try to e.g. load some more
    // libs and retry, or whatever seems relevant to react
    // to the error.
    
    return this->OnUnknownType(IdentifierType(id));
  }

  /// Print the list of creators we have.
  void print(std::ostream& os = std::cout) const 
  {
    typename CreatorMap::const_iterator it;
    for ( it = fCreators.begin(); it != fCreators.end(); ++it ) 
    {
      os << "Creator id=" << it->first 
	 << " for product " << it->second.productname()
	 << std::endl;
    }
  }

  /// Register a creator function for id.
  bool registerCreator(const IdentifierType& id, 
		       ProductCreator creator,
		       const char* productname)
  {
    bool ok = fCreators.insert(typename CreatorMap::value_type(id,ProductCreatorPair(creator,productname))).second;
    
    if (!ok)
    {
      std::cerr << "PHGenericFactoryT::registerCreator : registry of creator "
		<< "id " << id << " for product " << productname 
		<< "failed!" << std::endl;
    }
    else
    {
#ifdef DEBUG
      std::cout << "PHGenericFactoryT::registerCreator : creator id "
      << id << " for product " << productname 
      << " registered." << std::endl;
#endif
    }
    return ok;
  }

  /// Unregister a creator.
  bool unregisterCreator(const IdentifierType& id)
  {
    return fCreators.erase(id) == 1;
  }

private:

  PHGenericFactoryT() {}
  PHGenericFactoryT(const PHGenericFactoryT<AbstractProduct>&);
  PHGenericFactoryT<AbstractProduct>& operator=
  (const PHGenericFactoryT<AbstractProduct>&);
  ~PHGenericFactoryT() {}

  typedef std::map<IdentifierType,ProductCreatorPair> CreatorMap;
  CreatorMap fCreators;
};

#endif
