#include <geogram/basic/common.h>

#include <Mesh.hpp>







GEO::Mesh use_fTetWild()
{
    GEO::Initialize();

    Mesh mesh;
    Parameters& params = mesh.params;

    const unsigned int max_threads = std::numeric_limits<unsigned int>::max();
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads = std::min(max_threads, num_threads);
    params.num_threads = num_threads;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
}




