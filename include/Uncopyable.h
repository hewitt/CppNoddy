#ifndef UNCOPYABLE_H
#define UNCOPYABLE_H

namespace CppNoddy {

  /// An object to block copying

  class Uncopyable {

   protected:
    Uncopyable() {}
    ~Uncopyable() {}

   private:

    Uncopyable(const Uncopyable&);
    Uncopyable& operator=(const Uncopyable&);

  };

}

#endif // UNCOPYABLE_H
