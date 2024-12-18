import qa_functions
from qaoa_functions import *
import dask
from distributed import Scheduler,Worker
import asyncio
import argparse
import multiprocessing


def main():
    multiprocessing.set_start_method('spawn', force=True)

    parser = argparse.ArgumentParser(description="Input the size of the cluster.")
    parser.add_argument(
        "-n",
        "--number",
        type=int,
        choices=range(1, 65),
        required=True,
        help="Number of cores of the cluster (An integer from 1 to 64)",
    )
    args = parser.parse_args()
    num_workers = args.number

    print(f'Selected {num_workers} workers')

    asyncio.run(run_cluster(num_workers=num_workers))

async def run_cluster(num_workers):
    # Start the scheduler
    scheduler = await Scheduler(
        scheduler_file="phylo_schedule.json",
        port=8786,
        host="localhost",
        dashboard_address="127.0.0.1:8787",
    )
    print(f"Scheduler information saved to: phylo_schedule.json")

    # Start multiple workers
    workers = [await Worker(scheduler.address, nthreads=1) for _ in range(num_workers)]
    print(f"Workers connected: {[w.address for w in workers]}")

    # Keep the scheduler running
    await scheduler.finished()

if __name__ == "__main__":
    main()