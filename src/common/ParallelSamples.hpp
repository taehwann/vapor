#pragma once
#include <algorithm>
#include <future>
#include <thread>
#include <vector>

// Independent CPU samples only. Fixed chunks preserve per-sample arithmetic;
// futures propagate worker exceptions before the caller's data goes out of
// scope.
template <class Function> void parallelSamples(int count, Function &&sample) {
  int workers = std::min(4u, std::max(1u, std::thread::hardware_concurrency()));
  if (count < 1024 || workers == 1) {
    for (int i = 0; i < count; ++i)
      sample(i);
    return;
  }
  std::vector<std::future<void>> pending;
  auto chunk = [&](int worker) {
    for (int i = count * worker / workers; i < count * (worker + 1) / workers;
         ++i)
      sample(i);
  };
  for (int worker = 1; worker < workers; ++worker)
    pending.push_back(
        std::async(std::launch::async, [&, worker] { chunk(worker); }));
  chunk(0);
  for (auto &task : pending)
    task.get();
}
