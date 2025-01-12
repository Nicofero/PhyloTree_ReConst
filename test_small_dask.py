from dask.distributed import Client
from dask import delayed, compute
import time

def slow_function(x):
    time.sleep(1)
    return x

if __name__ == "__main__":
    client = Client()
    print(client)
    print(client.dashboard_link)
    time.sleep(10)
    tasks = [delayed(slow_function)(i) for i in range(5)]
    start_time = time.time()
    results = compute(*tasks)
    end_time = time.time()

    print(f"Results: {results}")
    print(f"Execution time: {end_time - start_time} seconds")
