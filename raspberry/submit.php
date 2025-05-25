<?php
$key = 'tajnykluczszyfr!';
$iv = '1234567890abcdef';
if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    $name = $_POST['name'] ?? 'brak';
    $surname = $_POST['surname'] ?? 'brak';
    $age = $_POST['age'] ?? 'brak';
    $rawMessage = $_POST['message'];
    $message = openssl_encrypt($rawMessage, 'AES-128-CBC', $key, 0, $iv);
    $time = date('Y-m-d H:i:s');

    $db = new SQLite3('responses.db');
    $db->exec('CREATE TABLE IF NOT EXISTS responses (
        id INTEGER PRIMARY KEY,
        name TEXT,
        surname TEXT,
        age TEXT,
        message TEXT,
        time TEXT
    )');

    $stmt = $db->prepare('INSERT INTO responses (name, surname, age, message, time) 
                          VALUES (:name, :surname, :age, :message, :time)');
    $stmt->bindValue(':name', $name, SQLITE3_TEXT);
    $stmt->bindValue(':surname', $surname, SQLITE3_TEXT);
    $stmt->bindValue(':age', $age, SQLITE3_TEXT);
    $stmt->bindValue(':message', $message, SQLITE3_TEXT);
    $stmt->bindValue(':time', $time, SQLITE3_TEXT);
    $stmt->execute();

    echo "Dane zapisane!";
} else {
    echo "Nieprawidłowe żądanie.";
}
