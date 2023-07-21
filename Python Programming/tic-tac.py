import logging
import random

logging.basicConfig(level=logging.INFO)

def print_board(board):
    """Prints the tic-tac-toe board with labeled columns.

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.
    """
    # Function to print the tic-tac-toe board with labeled columns
    print("\n    0   1   2 ")
    print("--------------")

    for i, row in enumerate(board):
        print(f"{i} | {row[0]} | {row[1]} | {row[2]}")
        if i < 2:
            print("--------------")

def check_winner(board, symbol):
    """Checks if a player with the given symbol has won the game.

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.
        symbol (str): The symbol ('X' or 'O') to check for a win.

    Returns:
        bool: True if the player with the given symbol has won, False otherwise.
    """
    # Function to check if a player with the given symbol has won the game
    for row in board:
        if all(cell == symbol for cell in row):  # Check rows
            return True

    for col in range(3):
        if all(board[row][col] == symbol for row in range(3)):  # Check columns
            return True

    if all(board[i][i] == symbol for i in range(3)):  # Check main diagonal
        return True

    if all(board[i][2 - i] == symbol for i in range(3)):  # Check secondary diagonal
        return True

    return False

def is_board_full(board):
    """Checks if the board is full (stalemate/tie).

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.

    Returns:
        bool: True if the board is full, False otherwise.
    """
    # Function to check if the board is full (stalemate/tie)
    return all(cell != " " for row in board for cell in row)

def get_empty_cells(board):
    """Gets a list of empty cells on the board.

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.

    Returns:
        list[tuple[int, int]]: A list of tuples representing the coordinates of empty cells.
    """
    # Function to get a list of empty cells on the board
    empty_cells = []
    for row in range(3):
        for col in range(3):
            if board[row][col] == " ":
                empty_cells.append((row, col))
    return empty_cells

def player_move(board, symbol):
    """Gets and validates the player's move.

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.
        symbol (str): The symbol ('X' or 'O') for the player.

    Raises:
        ValueError: If the input for row or column is not an integer.

    Returns:
        None
    """
    # Function to get and validate the player's move
    while True:
        try:
            row = int(input("Enter row (0, 1, or 2): "))
            col = int(input("Enter column (0, 1, or 2): "))
            if 0 <= row < 3 and 0 <= col < 3 and board[row][col] == " ":
                confirm = input(
                    f"Place '{symbol}' at row {row}, column {col}? (Yes/No): ").lower()
                if confirm == "yes":
                    board[row][col] = symbol
                    break
                else:
                    if(confirm)=="no".lower():
                        print("Move canceled.")

                    else:
                        print("Invalid Choice try again")
            else:
                logging.warning("Invalid move. Try again.")
        except ValueError:
            logging.warning("Invalid input. Try again.")

def computer_move(board, symbol):
    """Makes a random move for the computer.

    Args:
        board (list[list[str]]): The current state of the tic-tac-toe board.
        symbol (str): The symbol ('X' or 'O') for the computer.

    Returns:
        None
    """
    # Function for the computer to make a random move
    empty_cells = get_empty_cells(board)
    if empty_cells:
        row, col = random.choice(empty_cells)
        board[row][col] = symbol

def choose_player_symbol():
    """Asks the player for their preferred symbol (X or O).

    Returns:
        str: The chosen symbol ('X' or 'O').
    """
    # Function to ask the player for their preferred symbol (X or O)
    while True:
        player_symbol = input("Choose X or O: ").upper()
        if player_symbol in ["X", "O"]:
            return player_symbol
        else:
            print("Invalid choice. Please choose X or O.")

def play_game():
    """Plays a game of tic-tac-toe.

    Returns:
        None
    """
    # Main function to play the game
    board = [[" " for _ in range(3)] for _ in range(3)]
    player_symbol = choose_player_symbol()
    computer_symbol = "X" if player_symbol == "O" else "O"
    current_player = "X"  # Always start with 'X'

    print(
        f"You are '{player_symbol}' and the computer is '{computer_symbol}'.")
    print("Player X starts the game.")

    while not check_winner(board, player_symbol) and not check_winner(board, computer_symbol) and not is_board_full(board):
        print_board(board)
        print(f"Current player: {current_player}")

        if current_player == player_symbol:
            player_move(board, player_symbol)
        else:
            computer_move(board, computer_symbol)

        # Switch players
        current_player = player_symbol if current_player == computer_symbol else computer_symbol

    # Game over, print final board and result
    print_board(board)

    if check_winner(board, player_symbol):
        print(f"Congratulations! Player {player_symbol} wins!")
    elif check_winner(board, computer_symbol):
        print("Computer wins!")
    else:
        print("It's a tie!")

if __name__ == "__main__":
    play_game()
