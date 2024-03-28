#ifndef SYNC_H
#define SYNC_H

#include <atomic>

/* A busy waiting lock for accessing the history table and some other resources.
    - This is ideal over a mutex lock when access is expected to be a very short operation.
    - There may be room for OS / hardware specific preformance enchancements, but it seems unlikely to provide much benefit.
    - Does not prevent starvation but I doubt it will impact preformance. */
class spin_lock
{
public:
    /* Wait until the lock is released, then reserve it for yourself. */
    void lock()
    {
        while (lck.test_and_set(std::memory_order_acquire))
        {
            ;
        }
    };
    /* Release the lock for others to use. */
    void unlock()
    {
        lck.clear(std::memory_order_release);
    };

private:
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

#endif