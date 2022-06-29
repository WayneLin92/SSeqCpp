#ifndef THREAD_POOL_INCLUDED
#define THREAD_POOL_INCLUDED

#include <atomic>      // std::atomic
#include <functional>  // std::function
#include <memory>      // std::make_shared, std::make_unique, std::shared_ptr, std::unique_ptr
#include <mutex>       // std::mutex, std::scoped_lock, std::unique_lock
#include <queue>       // std::queue
#include <thread>      // std::thread

/**
 * The namespace `ut` provides basic utility functions
 */
namespace ut {

using concurrency_t = std::invoke_result_t<decltype(std::thread::hardware_concurrency)>;

class ThreadPool
{
private:
    std::vector<std::thread> threads_;
    std::queue<std::function<void()>> tasks_ = {};

    std::atomic<bool> running_ = false;
    std::atomic<bool> waiting_ = false;
    std::atomic<size_t> tasks_total_ = 0;
    std::mutex queue_mutex_ = {};
    std::mutex wait_mutex_ = {};
    std::condition_variable task_available_cv_ = {};
    std::condition_variable task_done_cv_ = {};

private:
    void worker()
    {
        while (running_) {
            std::function<void()> task;
            std::unique_lock<std::mutex> tasks_lock(queue_mutex_);
            task_available_cv_.wait(tasks_lock, [&] { return !tasks_.empty() || !running_; });
            if (running_) {
                task = std::move(tasks_.front());
                tasks_.pop();
                tasks_lock.unlock();
                task();
                --tasks_total_;
                if (waiting_) {
                    std::scoped_lock wait_lock(wait_mutex_);
                    task_done_cv_.notify_one();
                }
            }
        }
    }

public:
    explicit ThreadPool(const concurrency_t thread_count = std::thread::hardware_concurrency())
    {
        running_ = true;
        threads_.resize(thread_count);
        for (concurrency_t i = 0; i < thread_count; ++i)
            threads_[i] = std::thread(&ThreadPool::worker, this);
    }

    ~ThreadPool()
    {
        wait_for_tasks();
        running_ = false;
        task_available_cv_.notify_all();
        for (concurrency_t i = 0; i < threads_.size(); ++i)
            threads_[i].join();
    }

    void wait_for_tasks()
    {
        waiting_ = true;
        std::unique_lock<std::mutex> wait_lock(wait_mutex_);
        task_done_cv_.wait(wait_lock, [this] {
            return tasks_total_ == 0;
        });
        waiting_ = false;
    }

    template <typename Fn>
    void push_task(const Fn& task)
    {
        const std::scoped_lock tasks_lock(queue_mutex_);
        tasks_.push(std::function<void()>(task));
        ++tasks_total_;
        task_available_cv_.notify_one();
    }

    size_t get_threads_count() const
    {
        return threads_.size();
    }
};

}  // namespace ut

#endif