#
# Description:
# It compares two distinct strategies:
# 1. Mean Reversion
# 2. Breakout
# It provides separate, detailed summaries and trade logs for the best of each type.

import pandas as pd
import numpy as np
from datetime import timedelta

def run_reversion_backtest(df, entry_threshold, exit_threshold, max_hold_duration, tick_value_euro, stop_loss_euros, trade_size_lots, log_trades=False):

    signal_col = f'signal_rev_{entry_threshold}'
    entry_col = f'trade_entry_rev_{entry_threshold}'
    df[signal_col] = (df['deviation'].abs() > entry_threshold)
    df[entry_col] = (df[signal_col]) & (~df[signal_col].shift(1).fillna(False))
    entry_indices = df[df[entry_col]].index.to_list()
    
    total_pnl, trades_executed = 0.0, 0
    trade_log, trade_durations = [], []
    current_trade_idx, last_exit_bar_index = 0, -1

    while current_trade_idx < len(entry_indices):
        entry_index = entry_indices[current_trade_idx]
        if entry_index <= last_exit_bar_index:
            current_trade_idx += 1
            continue

        trades_executed += 1
        entry_price = df['spread'].loc[entry_index]
        entry_timestamp = df['Timestamp'].loc[entry_index]
        deviation_at_entry = df['deviation'].loc[entry_index]
        
        trade_type = 'SHORT' if deviation_at_entry > 0 else 'LONG'
        max_hold_exit_time = entry_timestamp + max_hold_duration
        stop_loss_points = stop_loss_euros / tick_value_euro
        stop_loss_price = entry_price + stop_loss_points if trade_type == 'SHORT' else entry_price - stop_loss_points

        future_df = df.loc[entry_index + 1:]
        
        reversion_indices = future_df.index[future_df['deviation'].abs() <= exit_threshold]
        profit_idx = reversion_indices.min() if not reversion_indices.empty else None

        time_exit_indices = future_df.index[future_df['Timestamp'] > max_hold_exit_time]
        time_exit_idx = time_exit_indices.min() if not time_exit_indices.empty else None

        sl_indices = future_df.index[future_df['spread'] >= stop_loss_price if trade_type == 'SHORT' else future_df['spread'] <= stop_loss_price]
        sl_idx = sl_indices.min() if not sl_indices.empty else None

        possible_exits = {'Profit Target': profit_idx, 'Time-Out': time_exit_idx, 'Stop-Loss': sl_idx}
        valid_exits = {reason: idx for reason, idx in possible_exits.items() if pd.notna(idx)}

        if not valid_exits:
            exit_reason, exit_idx = "End of Data", len(df) - 1
        else:
            exit_reason, exit_idx = min(valid_exits.items(), key=lambda item: item[1])

        exit_price = df['spread'].loc[exit_idx]
        exit_timestamp = df['Timestamp'].loc[exit_idx]
        last_exit_bar_index = exit_idx
        
        pnl_points = (entry_price - exit_price) if trade_type == 'SHORT' else (exit_price - entry_price)
        pnl_euros = pnl_points * tick_value_euro * trade_size_lots
        total_pnl += pnl_euros
        
        if log_trades:
            trade_duration = exit_timestamp - entry_timestamp
            trade_durations.append(trade_duration)
            trade_log.append({"EntryTimestamp": entry_timestamp, "EntryPrice": entry_price, "ExitTimestamp": exit_timestamp, "ExitPrice": exit_price, "ExitReason": exit_reason, "PnL_Euros": pnl_euros})

        current_trade_idx += 1
        
    if log_trades:
        return total_pnl, trades_executed, trade_log, trade_durations
    return total_pnl, trades_executed

def run_breakout_backtest(df, entry_threshold, breakout_tp_points, max_hold_duration, tick_value_euro, stop_loss_euros, trade_size_lots, log_trades=False):

    signal_col = f'signal_brk_{entry_threshold}'
    entry_col = f'trade_entry_brk_{entry_threshold}'
    df[signal_col] = (df['deviation'].abs() > entry_threshold)
    df[entry_col] = (df[signal_col]) & (~df[signal_col].shift(1).fillna(False))
    entry_indices = df[df[entry_col]].index.to_list()

    total_pnl, trades_executed = 0.0, 0
    trade_log, trade_durations = [], []
    current_trade_idx, last_exit_bar_index = 0, -1

    while current_trade_idx < len(entry_indices):
        entry_index = entry_indices[current_trade_idx]
        if entry_index <= last_exit_bar_index:
            current_trade_idx += 1
            continue

        trades_executed += 1
        entry_price = df['spread'].loc[entry_index]
        entry_timestamp = df['Timestamp'].loc[entry_index]
        deviation_at_entry = df['deviation'].loc[entry_index]

        trade_type = 'LONG' if deviation_at_entry > 0 else 'SHORT'
        max_hold_exit_time = entry_timestamp + max_hold_duration
        
        stop_loss_points = stop_loss_euros / tick_value_euro
        stop_loss_price = entry_price - stop_loss_points if trade_type == 'LONG' else entry_price + stop_loss_points
        profit_price = entry_price + breakout_tp_points if trade_type == 'LONG' else entry_price - breakout_tp_points

        future_df = df.loc[entry_index + 1:]

        profit_indices = future_df.index[future_df['spread'] >= profit_price if trade_type == 'LONG' else future_df['spread'] <= profit_price]
        profit_idx = profit_indices.min() if not profit_indices.empty else None

        sl_indices = future_df.index[future_df['spread'] <= stop_loss_price if trade_type == 'LONG' else future_df['spread'] >= stop_loss_price]
        sl_idx = sl_indices.min() if not sl_indices.empty else None
        
        time_exit_indices = future_df.index[future_df['Timestamp'] > max_hold_exit_time]
        time_exit_idx = time_exit_indices.min() if not time_exit_indices.empty else None

        possible_exits = {'Profit Target': profit_idx, 'Stop-Loss': sl_idx, 'Time-Out': time_exit_idx}
        valid_exits = {reason: idx for reason, idx in possible_exits.items() if pd.notna(idx)}

        if not valid_exits:
            exit_reason, exit_idx = "End of Data", len(df) - 1
        else:
            exit_reason, exit_idx = min(valid_exits.items(), key=lambda item: item[1])

        exit_price = df['spread'].loc[exit_idx]
        exit_timestamp = df['Timestamp'].loc[exit_idx]
        last_exit_bar_index = exit_idx

        pnl_points = (entry_price - exit_price) if trade_type == 'SHORT' else (exit_price - entry_price)
        pnl_euros = pnl_points * tick_value_euro * trade_size_lots
        total_pnl += pnl_euros
        
        if log_trades:
            trade_duration = exit_timestamp - entry_timestamp
            trade_durations.append(trade_duration)
            trade_log.append({"EntryTimestamp": entry_timestamp, "EntryPrice": entry_price, "ExitTimestamp": exit_timestamp, "ExitPrice": exit_price, "ExitReason": exit_reason, "PnL_Euros": pnl_euros})

        current_trade_idx += 1

    if log_trades:
        return total_pnl, trades_executed, trade_log, trade_durations
    return total_pnl, trades_executed

def main():
    # --- CONFIGURATION ---
    CSV_PATH = "strategy_data.csv"
    TICK_VALUE_EURO = 5.0
    TRADE_SIZE_LOTS = 1.0
    MAX_HOLD_HOURS = 8.0
    STOP_LOSS_EUROS = 100.0
    
    ENTRY_THRESHOLDS = np.arange(20, 101, 10)
    REVERSION_EXIT_THRESHOLDS = np.arange(5, 41, 5)
    BREAKOUT_TAKE_PROFIT_POINTS = np.arange(20, 101, 20)
    # --------------------------------------------------

    try:
        df = pd.read_csv(CSV_PATH, usecols=[0, 1, 2], header=0)
        df.columns = ['Timestamp', 'spread', 'ma']
        df['Timestamp'] = pd.to_datetime(df['Timestamp'], format='%d/%m/%Y %H:%M', errors='coerce')
        df.dropna(subset=['Timestamp'], inplace=True)
        df['spread'] = pd.to_numeric(df['spread'], errors='coerce')
        df['ma'] = pd.to_numeric(df['ma'], errors='coerce')
        df.sort_values(by='Timestamp', inplace=True)
        df.reset_index(drop=True, inplace=True)
        df['deviation'] = df['spread'] - df['ma']
    except Exception as e:
        print(f"Error loading data: {e}")
        return

    all_results = []
    hold_duration = timedelta(hours=MAX_HOLD_HOURS)

    print("\n--- TESTING REVERSION STRATEGIES (Trusted Logic) ---")
    for entry_thresh in ENTRY_THRESHOLDS:
        for rev_exit_thresh in REVERSION_EXIT_THRESHOLDS:
            pnl, trades = run_reversion_backtest(df.copy(), entry_thresh, rev_exit_thresh, hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS)
            if trades > 0:
                all_results.append({"type": "Reversion", "entry": entry_thresh, "exit": rev_exit_thresh, "exit_params_str": f"Exit Zone <= {rev_exit_thresh}", "total_pnl": pnl, "trades": trades, "avg_pnl": pnl / trades})
                print(f"Entry: {entry_thresh:3d} | Exit Zone <= {rev_exit_thresh:2d} | Trades: {trades:4d} | P&L: {pnl:10.2f} EUR")

    print("\n--- TESTING BREAKOUT STRATEGIES (Corrected Logic) ---")
    for entry_thresh in ENTRY_THRESHOLDS:
        for breakout_tp in BREAKOUT_TAKE_PROFIT_POINTS:
            pnl, trades = run_breakout_backtest(df.copy(), entry_thresh, breakout_tp, hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS)
            if trades > 0:
                all_results.append({"type": "Breakout", "entry": entry_thresh, "exit": breakout_tp, "exit_params_str": f"TP: {breakout_tp} pts", "total_pnl": pnl, "trades": trades, "avg_pnl": pnl / trades})
                print(f"Entry: {entry_thresh:3d} | TP: {breakout_tp:3d} pts | SL: {STOP_LOSS_EUROS:.0f} EUR | Trades: {trades:4d} | P&L: {pnl:10.2f} EUR")

    if not all_results:
        print("\nNo strategies produced trades.")
        return

    # --- Find Best of Each Type ---
    reversion_results = [res for res in all_results if res['type'] == 'Reversion']
    breakout_results = [res for res in all_results if res['type'] == 'Breakout']
    best_reversion = max(reversion_results, key=lambda x: x['total_pnl']) if reversion_results else None
    best_breakout = max(breakout_results, key=lambda x: x['total_pnl']) if breakout_results else None
    
    # --- Display Final Summaries ---
    if best_reversion:
        print("\n\n--- Best Performing REVERSION Strategy ---")
        print(f"  - Entry Threshold: {best_reversion['entry']} points")
        print(f"  - Exit Rules:      {best_reversion['exit_params_str']}")
        print(f"  - Total P&L:       {best_reversion['total_pnl']:.2f} EUR")
        print(f"  - Average P&L:     {best_reversion['avg_pnl']:.2f} EUR per trade")
        _, _, _, durations = run_reversion_backtest(df.copy(), best_reversion['entry'], best_reversion['exit'], hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS, log_trades=True)
        if durations:
            avg_duration = pd.to_timedelta(durations).mean()
            print(f"  - Avg Time in Trade: {str(avg_duration).split('.')[0]}")
    
    if best_breakout:
        print("\n--- Best Performing BREAKOUT Strategy ---")
        print(f"  - Entry Threshold: {best_breakout['entry']} points")
        print(f"  - Exit Rules:      {best_breakout['exit_params_str']}")
        print(f"  - Total P&L:       {best_breakout['total_pnl']:.2f} EUR")
        print(f"  - Average P&L:     {best_breakout['avg_pnl']:.2f} EUR per trade")
        _, _, _, durations = run_breakout_backtest(df.copy(), best_breakout['entry'], best_breakout['exit'], hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS, log_trades=True)
        if durations:
            avg_duration = pd.to_timedelta(durations).mean()
            print(f"  - Avg Time in Trade: {str(avg_duration).split('.')[0]}")

    # --- Log Trades for Both Winners ---
    print("\nGenerating trade logs for the best strategies...")
    if best_reversion:
        _, _, reversion_log, _ = run_reversion_backtest(df.copy(), best_reversion['entry'], best_reversion['exit'], hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS, log_trades=True)
        pd.DataFrame(reversion_log).to_csv("reversion_trade_log.csv", index=False)
        print("  - Saved reversion_trade_log.csv")
    if best_breakout:
        _, _, breakout_log, _ = run_breakout_backtest(df.copy(), best_breakout['entry'], best_breakout['exit'], hold_duration, TICK_VALUE_EURO, STOP_LOSS_EUROS, TRADE_SIZE_LOTS, log_trades=True)
        pd.DataFrame(breakout_log).to_csv("breakout_trade_log.csv", index=False)
        print("  - Saved breakout_trade_log.csv")

if __name__ == "__main__":
    main()

