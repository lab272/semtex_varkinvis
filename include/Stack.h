/*****************************************************************************
 * Stack.h: operations for LIFO stack.                                       *
 *                                                                           *
 * Summary:                                                                  *
 *   creator, destructor, push, pop, isEmpty.                                *
 *****************************************************************************/

/* $Id$ */

#ifndef StackH
#define StackH

#include <iostream.h>


template<class T>
class Stack {
public:
  Stack  () : top (0) {}
  ~Stack ();

  void  push    (T x);
  T     pop     ();
  int   isEmpty ();

private:
  class Node {
  public:
    Node (T x) : next (0), datum (x) {}
    Node* next;
    T     datum;
  };

  Node* top;

  Stack(const Stack<T>&);                // Prohibit, since not implemented.
  Stack<T>& operator=(const Stack<T>&);  // Prohibit, since not implemented.
};





template<class T>
inline void Stack<T>::push(T x) {
  if (isEmpty())
    top = new Node (x);
  else {
    Node *p = new Node (x);
    p -> next = top;
    top = p;
  }
}


template<class T>
inline T Stack<T>::pop() {
  if (isEmpty()) {
    cout << "WARNING: stack empty" << endl;
    return 0;
  } else {
    Node *p = top;
    T     value = top -> datum;
    top = top -> next;
    delete p;
    return value;
  }
}


template<class T>
inline int Stack<T>::isEmpty() {
  return top == 0;
}


template<class T>
inline Stack<T>::~Stack() {
  while (! isEmpty()) {
    Node* p = top -> next;
    delete top;
    top = p;
  }
}


#endif