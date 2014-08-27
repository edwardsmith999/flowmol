/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_SERIAL_TASKPOLICY_HPP
#define KOKKOS_SERIAL_TASKPOLICY_HPP

#include <string>
#include <typeinfo>
#include <stdexcept>

#include <Kokkos_Serial.hpp>
#include <Kokkos_TaskPolicy.hpp>
#include <Kokkos_View.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class TaskMember< Kokkos::Serial , void , void >
{
public:

  friend class TaskManager< Kokkos::Serial > ;

  enum { MAX_DEPENDENCE = 13 };

  /**\brief  States of a task */
  enum { STATE_CONSTRUCTING = 0 , STATE_WAITING = 1 , STATE_EXECUTING = 2 , STATE_COMPLETE = 4 };

  /**\brief  Base dependence count when a task is allocated.
   *         A separate dependence array is allocated when the number
   *         of dependences exceeds this count.
   */

  typedef void (* function_type)( TaskMember * );

  const std::type_info & m_typeid ;
  const function_type    m_dealloc ;
  const function_type    m_apply ;

private:

  int            m_state ;
  int            m_ref_count ; ///< Reference count
  TaskMember   * m_wait ;      ///< Linked list of tasks waiting on this task.
  TaskMember   * m_next ;      ///< This task is a member of a linked list of
                               ///< tasks waiting on another task.
  TaskMember   * m_dep[ MAX_DEPENDENCE ]; ///< Dependences of this task

  TaskMember( const TaskMember & );
  TaskMember & operator = ( const TaskMember & );

protected :

  inline
  TaskMember( const function_type    arg_dealloc
            , const function_type    arg_apply
            , const std::type_info & arg_type = typeid(void)
            )
    : m_typeid(  arg_type )
    , m_dealloc( arg_dealloc )
    , m_apply(   arg_apply )
    , m_state( STATE_CONSTRUCTING )
    , m_ref_count(0)
    , m_wait(0)
    , m_next(0)
    { for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) m_dep[i] = 0 ; }

public:

  template < class DerivedTaskMember >
  static
  void deallocate( TaskMember * t )
    { delete static_cast< DerivedTaskMember * >( t ); }

  inline static
  TaskMember * verify_type( TaskMember * t ) { return t ; }

  typedef FutureValueTypeIsVoidError get_result_type ;

  get_result_type get() const { return get_result_type() ; }

  inline
  TaskMember * get_dependence( int i ) const
    { return ( STATE_EXECUTING == m_state && 0 <= i && i < MAX_DEPENDENCE ) ? m_dep[i] : (TaskMember*) 0 ; }

  inline
  int get_dependence() const
    {
      int i = 0 ;
      if ( STATE_EXECUTING == m_state ) { for ( ; i < MAX_DEPENDENCE && m_dep[i] != 0 ; ++i ); }
      return i ;
    }
};

//----------------------------------------------------------------------------

template<>
class TaskManager< Kokkos::Serial >
{
public:

  typedef TaskMember< Kokkos::Serial > task_root_type ;

  enum { MAX_DEPENDENCE = task_root_type::MAX_DEPENDENCE };

  void schedule( task_root_type * );

  static void verify_set_dependence( task_root_type * , int );

  static void assign( task_root_type ** const , task_root_type * );

  void wait( task_root_type * );

  TaskManager();
  TaskManager( const TaskManager & );
  TaskManager & operator = ( const TaskManager & );


  template< class A1 , class A2 >
  void wait( const Future<A1,A2> & future ) { wait( future.m_task ); }

  template< class A1 , class A2 >
  void set_dependence( task_root_type * t
                     , const Future<A1,A2> * const dep
                     , typename Impl::enable_if
                        < Impl::is_same< typename Future<A1,A2>::execution_space , Kokkos::Serial >::value
                        , const int >::type n
                      )
    {
      verify_set_dependence( t , n );
      int i = 0 ;
      for ( ; i < n ; ++i ) assign( & t->m_dep[i] , dep[i].m_task );
      for ( ; i < MAX_DEPENDENCE ; ++i ) assign( & t->m_dep[i] , 0 );
    }

private:

  task_root_type * m_ready ;
  task_root_type * const m_denied ;

};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template < class ResultType >
class TaskMember< Kokkos::Serial , ResultType , void > : public TaskMember< Kokkos::Serial >
{
private:

protected:

  typedef TaskMember< Kokkos::Serial >::function_type  function_type ;

  inline
  TaskMember( const function_type    arg_dealloc
            , const function_type    arg_apply
            )
    : TaskMember< Kokkos::Serial >( arg_dealloc , arg_apply , typeid(ResultType) )
    , m_result()
    {}

public:

  ResultType  m_result ;

  inline static
  TaskMember *
  verify_type( TaskMember< Kokkos::Serial > * t )
    {
      if ( t != 0 && t->m_typeid != typeid(ResultType) ) {
        throw std::runtime_error( std::string("Kokkos::Future bad cast for result type"));
      }
      return static_cast< TaskMember *>( t );
    }

  typedef const ResultType & get_result_type ;

  inline
  get_result_type get() const { return m_result ; }
};

//----------------------------------------------------------------------------

template< class ResultType , class FunctorType >
class TaskMember< Kokkos::Serial , ResultType , FunctorType >
  : public TaskMember< Kokkos::Serial , ResultType >
  , public FunctorType
{
private:

  typedef TaskMember< Kokkos::Serial >              member_root_type ;
  typedef TaskMember< Kokkos::Serial , ResultType > member_base_type ;

  static
  void apply( member_root_type * t )
    {
      member_base_type * m = static_cast< member_base_type * >(t);
      static_cast< TaskMember * >(m)->FunctorType::apply( m->m_result );
    }

protected:

  inline 
  TaskMember( const typename member_root_type::function_type  arg_dealloc
            , const typename member_root_type::function_type  arg_apply
            , const FunctorType &  arg_functor
            )
    : member_base_type( arg_dealloc , arg_apply )
    , FunctorType( arg_functor )
    {}

public:

  inline 
  TaskMember( const FunctorType &  arg_functor )
    : member_base_type( & member_root_type::template deallocate< TaskMember >
                      , & TaskMember::apply )
    , FunctorType( arg_functor )
    {}
};

//----------------------------------------------------------------------------

template< class FunctorType >
class TaskMember< Kokkos::Serial , void , FunctorType >
  : public TaskMember< Kokkos::Serial >
  , public FunctorType
{
private:

  typedef TaskMember< Kokkos::Serial > member_root_type ;

  static
  void apply( member_root_type * t )
    { static_cast< TaskMember * >(t)->FunctorType::apply(); }

protected:

  inline 
  TaskMember( const typename member_root_type::function_type  arg_dealloc
            , const typename member_root_type::function_type  arg_apply
            , const FunctorType &  arg_functor
            )
    : member_root_type( arg_dealloc , arg_apply )
    , FunctorType( arg_functor )
    {}

public:

  inline 
  TaskMember( const FunctorType &  arg_functor )
    : member_root_type( & member_root_type::template deallocate< TaskMember >
                      , & TaskMember::apply )
    , FunctorType( arg_functor )
    {}
};

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

//----------------------------------------------------------------------------

template<>
class TaskPolicy< Impl::TaskDepends< Kokkos::Serial > >
{
public:

  typedef Kokkos::Serial execution_space ;

private:

  enum { MAX_DEPENDENCE = Impl::TaskMember< execution_space >::MAX_DEPENDENCE };

  Kokkos::Impl::TaskManager< execution_space >  & m_task_manager ;
  Kokkos::Future< execution_space >               m_depends[ MAX_DEPENDENCE ];

  TaskPolicy();
  TaskPolicy & operator = ( const TaskPolicy & );

public:

  template< typename A1 , typename A2 >
  TaskPolicy( Kokkos::Impl::TaskManager< execution_space > & manager
            , const size_t n
            , const Future< A1 , A2 > * const dep )
    : m_task_manager( manager )
    {
      int i = 0 ;
      for ( ; i < n ; ++i ) m_depends[i] = dep[i] ;
      for ( ; i < MAX_DEPENDENCE ; ++i ) m_depends[i] = Future< execution_space >();
    }

  // Spawn a serial task:
  template< class FunctorType , class ValueType >
  Future< ValueType , execution_space >
  spawn( const FunctorType & functor ) const
    {
      // Allocate a copy functor and insert into queue
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType > member_type ;
      member_type * m = new member_type( functor );
      m_task_manager.set_dependence( m , m_depends );
      m_task_manager.schedule( m );
      return Future< ValueType , execution_space >( m );
    }

  // Construct a task policy for foreach-range tasks:
  // spawn( task_policy.depends(N,d).foreach(RangePolicy) , functor );
  // spawn( task_policy.foreach(RangePolicy) , functor );
  template< class ExecPolicy >
  TaskPolicy< Impl::TaskForEach< ExecPolicy > >
  foreach( const ExecPolicy & arg_policy )
    { return TaskPolicy< Impl::TaskForEach< ExecPolicy > >( m_task_manager , arg_policy , m_depends ); }

  // Construct a task policy for reduce-range tasks:
  template< class ExecPolicy >
  TaskPolicy< Impl::TaskForEach< ExecPolicy > >
  reduce( const ExecPolicy & arg_policy )
    { return TaskPolicy< Impl::TaskReduce< ExecPolicy > >( m_task_manager , arg_policy , m_depends ); }
};

//----------------------------------------------------------------------------

template<>
class TaskPolicy< Kokkos::Serial >
{
public:

  typedef Kokkos::Serial execution_space ;

private:

  typedef Impl::TaskMember< execution_space , void , void > task_base_type ;

  Kokkos::Impl::TaskManager< execution_space > & m_task_manager ;

  template< class FunctorType >
  static
  void apply( task_base_type * t )
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType >  member_type ;
      static_cast< member_type * >(t)->FunctorType::apply();
    }

  TaskPolicy & operator = ( const TaskPolicy & );

public:

  TaskPolicy();
  TaskPolicy( const TaskPolicy & rhs )
    : m_task_manager( rhs.m_task_manager ) {}

  // Requires:
  // class DerivedMemberType : public TaskMember< execution_space , typename FunctorType::value_type , FunctorType > ...
  template< class FunctorType >
  Future< void , execution_space >
  get_dependence( const FunctorType * task_functor , int i ) const
    {
      typedef const Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType >  member_type ;
      return Future<void,execution_space>( static_cast< member_type * >(task_functor)->task_base_type::get_dependence(i) );
    }

  template< class FunctorType >
  int get_dependence( const FunctorType * task_functor ) const
    {
      typedef const Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType >  member_type ;
      return static_cast< member_type * >(task_functor)->task_base_type::get_dependence();
    }

  template< class A1 , class A2 >
  void wait( const Future<A1,A2> & future ) const { m_task_manager.wait( future ); }

  template< class FunctorType , class A1 , class A2 >
  void respawn( FunctorType * task_functor
              , const Future<A1,A2> * const dep
              , typename Impl::enable_if
                  < Impl::is_same< typename Future<A1,A2>::execution_space , execution_space >::value
                  , const int
                  >::type n
              ) const
    {
      typedef Impl::TaskMember< execution_space , typename FunctorType::value_type , FunctorType >  member_type ;
      m_task_manager.set_dependence( static_cast< member_type * >( task_functor ) , dep , n );
      m_task_manager.schedule( static_cast< member_type * >( task_functor ) );
    }

  // Allocate a copy functor and insert into queue
  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  spawn( const FunctorType & functor ) const
    {
      typedef typename FunctorType::value_type value_type ;
      typedef Impl::TaskMember< execution_space , value_type , FunctorType >  member_type ;
      member_type * m = new member_type( functor );
      m_task_manager.schedule( m );
      return Future< value_type , execution_space >( m );
    }

  // Construct a task policy with dependences:
  // spawn( task_policy.depends(N,d) , functor );
  template< class A1 , class A2 >
  TaskPolicy< Impl::TaskDepends< execution_space > >
  depends( const Future< A1 , A2 > * const d
         , typename Impl::enable_if<
             ( Impl::is_same< typename Future<A1,A2>::execution_space , execution_space >::value
             ), const int >::type n 
         )
    { return TaskPolicy< Impl::TaskDepends< execution_space > >( m_task_manager , n , d ); }

  // Construct a task policy for foreach-range tasks:
  // spawn( task_policy.depends(N,d).foreach(RangePolicy) , functor );
  // spawn( task_policy.foreach(RangePolicy) , functor );
  template< class ExecPolicy >
  TaskPolicy< Impl::TaskForEach< ExecPolicy > >
  foreach( const ExecPolicy & arg_policy )
    { return TaskPolicy< Impl::TaskForEach< ExecPolicy > >( m_task_manager , arg_policy ); }

  // Construct a task policy for reduce-range tasks:
  template< class ExecPolicy >
  TaskPolicy< Impl::TaskReduce< ExecPolicy > >
  reduce( const ExecPolicy & arg_policy )
    { return TaskPolicy< Impl::TaskReduce< ExecPolicy > >( m_task_manager , arg_policy ); }
};

//----------------------------------------------------------------------------

template< typename IntType , unsigned P >
class TaskPolicy< Impl::TaskForEach< Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P >  >  >
{
public:

  typedef Kokkos::Serial execution_space ;

private:

  typedef RangePolicy< execution_space , void , IntType , P > range_policy ;
  typedef Impl::TaskManager< execution_space >  task_manager ;
  typedef Impl::TaskMember<  execution_space >  task_root_type ;

  task_manager & m_task_manager ;
  range_policy   m_range_policy ;

  // ForEach task
  template< class FunctorType >
  class member_type : public Impl::TaskMember< Kokkos::Serial , void , FunctorType >
  {
  private:

    typedef Impl::TaskMember< Kokkos::Serial , void , FunctorType >    task_base_type ;

    range_policy  m_policy ;

    static
    void apply( task_root_type * t )
      {
        range_policy const & r  = static_cast< member_type * >( static_cast< task_base_type * >( t ) )->m_policy ;
        FunctorType        & f  = * static_cast< FunctorType * >( static_cast< task_base_type * >( t ) );
        FunctorType  const & cf = f ;

        const IntType e = r.end();
        for ( IntType i = r.begin() ; i < e ; ++i ) { cf(i); }
        f.apply();
      }

  public:

    member_type( const FunctorType  & arg_func 
               , const range_policy & arg_policy
               )
      : task_base_type( & task_base_type::template deallocate< member_type >
                      , & member_type::apply
                      , arg_func
                      )
      , m_policy( arg_policy )
      {}
  };


  TaskPolicy();
  TaskPolicy & operator = ( const TaskPolicy & );

public:

  TaskPolicy( task_manager & manager , const range_policy & policy )
    : m_task_manager( manager )
    , m_range_policy( policy )
    {}

  template< class FunctorType , class ValueType >
  Future< ValueType , execution_space >
  spawn( const FunctorType & functor ) const
    {
      typedef Future< ValueType , execution_space > future_type ;

      // Allocate a copy functor and insert into queue

      task_root_type * const t = new member_type< FunctorType >( functor , m_range_policy );

      m_task_manager.schedule( t );

      return future_type( t );
    }
};

//----------------------------------------------------------------------------

template< typename IntType , unsigned P >
class TaskPolicy< Impl::TaskReduce< Kokkos::RangePolicy< Kokkos::Serial , void , IntType , P >  >  >
{
public:

  typedef Kokkos::Serial execution_space ;

private:

  typedef RangePolicy< execution_space , void , IntType , P >  range_policy ;
  typedef Impl::TaskManager< execution_space >  task_manager ;
  typedef Impl::TaskMember<  execution_space >  task_root_type ;

  task_manager & m_task_manager ;
  range_policy   m_range_policy ;

  // ForEach task
  template< class FunctorType >
  class member_type : public Impl::TaskMember< Kokkos::Serial , typename FunctorType::value_type , FunctorType >
  {
  private:
    typedef typename FunctorType::value_type value_type ;

    typedef Impl::TaskMember< Kokkos::Serial , value_type , FunctorType >    task_base_type ;
    typedef Impl::TaskMember< Kokkos::Serial , value_type >    task_value_type ;

    range_policy  m_policy ;

    static
    void apply( task_root_type * t )
      {
        task_base_type     & b  = * static_cast< task_base_type * >( t );
        range_policy const & r  = static_cast< member_type & >( b ).m_policy ;
        FunctorType        & f  = static_cast< FunctorType & >( b );
        FunctorType  const & cf = f ;

        cf.init( b.m_result );
        const IntType e = r.end();
        for ( IntType i = r.begin() ; i < e ; ++i ) { cf(i,b.m_result); }
        f.apply( b.m_result );
      }

  public:

    member_type( const FunctorType  & arg_func 
               , const range_policy & arg_policy
               )
      : task_base_type( & task_base_type::template deallocate< member_type >
                      , & member_type::apply
                      , arg_func
                      )
      , m_policy( arg_policy )
      {}
  };

  TaskPolicy();
  TaskPolicy & operator = ( const TaskPolicy & );

public:

  TaskPolicy( task_manager & manager , const range_policy & policy )
    : m_task_manager( manager )
    , m_range_policy( policy )
    {}

  template< class FunctorType >
  Future< typename FunctorType::value_type , execution_space >
  spawn( const FunctorType & functor ) const
    {
      typedef Future< typename FunctorType::value_type , execution_space > future_type ;

      // Allocate a copy functor and insert into queue

      task_root_type * const t = new member_type<FunctorType>( functor , m_range_policy );

      m_task_manager.schedule( t );

      return future_type( t );
    }
};

} // namespace Kokkos

//----------------------------------------------------------------------------

#endif /* #define KOKKOS_SERIAL_TASK_HPP */

