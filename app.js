const DATA_URL = "data/words.json";

let allCards = [];
let cards = [];
let baseSessionCards = [];
let sessionAnswers = [];
let currentIndex = 0;
let correctCount = 0;
let wrongCount = 0;

let answerVisible = false;
let isChangingCard = false;
let writeComparisonDone = false;
let stats = {};

let isNewWordsMode = false;
let newWordsPool = [];
let newWordsActiveCount = 1;
let newWordsRoundResults = {};
let allNewWordsIntroduced = false;

let heldNewWordIds = new Set();
let isManagingNewWordsDeck = false;

// Controle especial das palavras acertadas na primeira tentativa
let newWordsCardProgress = {};

// Quantas rodadas uma palavra acertada de primeira ficará afastada
const NEW_WORD_FIRST_TRY_DELAY = 15;

let sessionStartTime = null;
let elapsedBeforePause = 0;
let timerIntervalId = null;

let wordsDirection = "sv-pt";
let wordsViewMode = "all";
let expandedWordCardId = null;
let expandedWordsLetter = null;
let cardsSinceUnseen = 0;

let currentExercise = null;
let exerciseFinished = false;
let exerciseErrorsCopyText = "";

let exerciseVocabularyIndex = new Map();
let expandedExerciseVocabulary = null;

let exerciseConsultedVocabulary = new Map();

const STORAGE_KEY = "flashcardsSuecoStats";
const THEME_STORAGE_KEY = "flashcardsSuecoTheme";
const HELD_NEW_WORDS_STORAGE_KEY = "flashcardsSuecoHeldNewWords";

const UNDEREXPOSURE_PRIORITY_PER_VIEW = 25;
const UNDEREXPOSURE_MAX_GAP = 6;
const UNSEEN_CARD_MAX_GAP = 4;

const UNSEEN_MIN_ACCURACY = 0.60;
const UNSEEN_ACCURACY_WINDOW = 10;
const UNSEEN_MIN_ANSWERS = 5;

const UNSEEN_HARD_MAX_GAP = 15;
const UNSEEN_HARD_MIN_ACCURACY = 0.50;

const PRONUNCIATION_RULES = [{
    title: "1. A vogal depois de SK, G e K muda o som",
    text: "Quando SK, G ou K vêm antes de A, O, U ou Å, o som tende a ser mais “duro”. Antes de E, I, Y, Ä ou Ö, o som costuma ficar mais “suave”.",
    examples: ["gott", "gift", "gör", "göra", "kul", "komma", "kort", "att känna", "att köpa", "skild"]
  },
  {
    title: "2. Combinações que não soam como seriam lidas em português",
    text: "Algumas combinações têm som especial: TJ soa como um chiado suave; SJ, SKJ, STJ e SK antes de vogal suave têm um som de “ch” mais forte; CH varia conforme a palavra.",
    examples: ["Tjena!", "sju", "att sjunga", "skild"]
  },
  {
    title: "3. RS, RT, RD e RN se juntam na fala",
    text: "Depois de R, as letras S, T, D e N costumam ser pronunciadas com a língua mais para trás. Na prática, não se fala como duas letras totalmente separadas.",
    examples: ["torsdag", "Ursäkta?", "ett barn", "snart", "kort", "ett hjärta"]
  },
  {
    title: "4. Algumas letras somem ou enfraquecem na fala rápida",
    text: "O G final em palavras terminadas em -dag muitas vezes quase desaparece. Em várias palavras frequentes, letras como R, G ou H podem ficar bem fracas dependendo da frase e da região.",
    examples: ["måndag", "tisdag", "onsdag", "fredag", "lördag", "söndag"]
  },
  {
    title: "5. A sílaba tônica normalmente fica na parte mais importante da palavra",
    text: "Palavras suecas comuns muitas vezes têm a primeira sílaba forte. Palavras com be- ou för- podem puxar a força para a segunda sílaba. Em frases, substantivos, verbos, adjetivos e advérbios costumam receber mais ênfase.",
    examples: ["en bok", "att äta", "förstå", "forska", "komma"]
  },
  {
    title: "6. Vogal longa x consoante longa",
    text: "Uma vogal seguida de uma só consoante costuma soar mais longa. Quando há duas ou mais consoantes depois, a vogal costuma ficar mais curta e a consoante seguinte pesa mais.",
    examples: ["en bok", "noll", "åtta", "komma", "kort"]
  },
  {
    title: "7. Melodia da frase",
    text: "O sueco não é falado “reto”: a voz sobe e desce, especialmente em vogais longas. Em afirmações, a melodia geralmente cai no final; em perguntas, pode subir ou variar conforme a região.",
    examples: ["Hej!", "Tjena!", "Ursäkta?", "Jag förstår inte."]
  }
];

const newWordsToggleButton = document.querySelector("#newWordsToggleButton");

const newWordsDeckActions =
  document.querySelector("#newWordsDeckActions");

const holdCurrentWordButton =
  document.querySelector("#holdCurrentWordButton");

const restoreHeldWordsButton =
  document.querySelector("#restoreHeldWordsButton");

const summaryScreen = document.querySelector("#summaryScreen");

const themeToggleButton = document.querySelector("#themeToggleButton");

const summaryTime = document.querySelector("#summaryTime");

const finishSessionButton = document.querySelector("#finishSessionButton");
const repeatSessionButton = document.querySelector("#repeatSessionButton");
const newSessionButton = document.querySelector("#newSessionButton");

const summaryTotal = document.querySelector("#summaryTotal");
const summaryCorrect = document.querySelector("#summaryCorrect");
const summaryWrong = document.querySelector("#summaryWrong");
const summaryCorrectRate = document.querySelector("#summaryCorrectRate");
const summaryWrongList = document.querySelector("#summaryWrongList");

const setupScreen = document.querySelector("#setupScreen");
const studyScreen = document.querySelector("#studyScreen");

const typeFilterGroup = document.querySelector("#typeFilterGroup");
const themeFilterGroup = document.querySelector("#themeFilterGroup");
const sourceFilterGroup = document.querySelector("#sourceFilterGroup");
const chapterFilterLabel = document.querySelector("#chapterFilterLabel");
const chapterFilterGroup = document.querySelector("#chapterFilterGroup");
const sourceTitleFilterLabel = document.querySelector("#sourceTitleFilterLabel");
const sourceTitleFilterGroup = document.querySelector("#sourceTitleFilterGroup");

const startSessionButton = document.querySelector("#startSessionButton");
const backToSetupButton = document.querySelector("#backToSetupButton");
const setupMessage = document.querySelector("#setupMessage");

const pronunciationRulesButton = document.querySelector("#pronunciationRulesButton");
const pronunciationScreen = document.querySelector("#pronunciationScreen");
const backFromPronunciationButton = document.querySelector("#backFromPronunciationButton");
const pronunciationRulesList = document.querySelector("#pronunciationRulesList");

const grammarButton = document.querySelector("#grammarButton");
const grammarScreen = document.querySelector("#grammarScreen");
const backFromGrammarButton = document.querySelector("#backFromGrammarButton");

const pluralRulesButton = document.querySelector("#pluralRulesButton");
const pluralScreen = document.querySelector("#pluralScreen");
const backFromPluralButton = document.querySelector("#backFromPluralButton");

const verbFormsButton = document.querySelector("#verbFormsButton");
const verbFormsScreen = document.querySelector("#verbFormsScreen");
const backFromVerbFormsButton = document.querySelector("#backFromVerbFormsButton");

const resetStatsButton = document.querySelector("#resetStatsButton");

const wordsButton = document.querySelector("#wordsButton");
const wrongWordsButton = document.querySelector("#wrongWordsButton");
const wordsScreen = document.querySelector("#wordsScreen");
const backFromWordsButton = document.querySelector("#backFromWordsButton");
const wordsScreenTitle = document.querySelector("#wordsScreenTitle");
const wordsScreenDescription = document.querySelector("#wordsScreenDescription");
const wordsDirectionButton = document.querySelector("#wordsDirectionButton");
const wordsCount = document.querySelector("#wordsCount");
const wordsList = document.querySelector("#wordsList");

const exerciseButton = document.querySelector("#exerciseButton");
const exerciseScreen = document.querySelector("#exerciseScreen");
const backFromExerciseButton = document.querySelector("#backFromExerciseButton");

const exerciseImportPanel = document.querySelector("#exerciseImportPanel");
const exerciseSourceInput = document.querySelector("#exerciseSourceInput");
const generateExerciseButton = document.querySelector("#generateExerciseButton");
const exerciseImportMessage = document.querySelector("#exerciseImportMessage");

const exerciseRenderedPanel = document.querySelector("#exerciseRenderedPanel");
const exerciseTitle = document.querySelector("#exerciseTitle");
const exerciseResultSummary = document.querySelector("#exerciseResultSummary");
const exerciseBlocks = document.querySelector("#exerciseBlocks");

const finishExerciseButton = document.querySelector("#finishExerciseButton");
const copyExerciseErrorsButton = document.querySelector("#copyExerciseErrorsButton");
const newExerciseButton = document.querySelector("#newExerciseButton");

const appScreens = [
  setupScreen,
  studyScreen,
  summaryScreen,
  pronunciationScreen,
  grammarScreen,
  pluralScreen,
  verbFormsScreen,
  wordsScreen,
  exerciseScreen
];

function hideAllScreens() {
  appScreens.forEach((screen) => {
    screen.classList.add("hidden");
  });
}

function showScreen(screen) {
  hideAllScreens();
  screen.classList.remove("hidden");

  newWordsToggleButton.classList.toggle(
    "hidden",
    screen !== setupScreen
  );

  const shouldShowDeckActions =
    screen === studyScreen &&
    isNewWordsMode;

  newWordsDeckActions.classList.toggle(
    "hidden",
    !shouldShowDeckActions
  );

  updateHeldNewWordsUI();
}

const questionLabel = document.querySelector("#questionLabel");

const directionSelect = document.querySelector("#directionSelect");
const answerModeSelect = document.querySelector("#answerModeSelect");
const answerModeLabel = document.querySelector("#answerModeLabel");

const answerLabel = document.querySelector("#answerLabel");
const backHint = document.querySelector("#backHint");

const writeBox = document.querySelector("#writeBox");
const answerInput = document.querySelector("#answerInput");
const checkAnswerButton = document.querySelector("#checkAnswerButton");
const nextWriteCardButton = document.querySelector("#nextWriteCardButton");
const writeResultBox = document.querySelector("#writeResultBox");
const userAnswerText = document.querySelector("#userAnswerText");
const comparisonFeedbackText = document.querySelector("#comparisonFeedbackText");

const cardAudioWrap = document.querySelector("#cardAudioWrap");
const cardAudio = document.querySelector("#cardAudio");
const playAudioButton = document.querySelector("#playAudioButton");

const showSwedishFromAudioButton = document.querySelector("#showSwedishFromAudioButton");
const audioSwedishText = document.querySelector("#audioSwedishText");

const questionText = document.querySelector("#questionText");
const flashcard = document.querySelector("#flashcard");
const answerText = document.querySelector("#answerText");
const answerGrammarForm = document.querySelector("#answerGrammarForm");

const pronunciationButton = document.querySelector("#pronunciationButton");

const frontPronunciationButton = document.querySelector("#frontPronunciationButton");
const backPronunciationButton = document.querySelector("#backPronunciationButton");

const resultButtons = document.querySelector("#resultButtons");
const correctButton = document.querySelector("#correctButton");
const wrongButton = document.querySelector("#wrongButton");
const message = document.querySelector("#message");

function loadHeldNewWords() {
  try {
    const saved = localStorage.getItem(
      HELD_NEW_WORDS_STORAGE_KEY
    );

    if (!saved) {
      heldNewWordIds = new Set();
      return;
    }

    const parsed = JSON.parse(saved);

    heldNewWordIds = new Set(
      Array.isArray(parsed)
        ? parsed.map(String)
        : []
    );
  } catch (error) {
    console.error(
      "Erro ao carregar palavras guardadas:",
      error
    );

    heldNewWordIds = new Set();

    try {
      localStorage.removeItem(
        HELD_NEW_WORDS_STORAGE_KEY
      );
    } catch (storageError) {
      console.error(
        "Erro ao limpar palavras guardadas:",
        storageError
      );
    }
  }
}

function saveHeldNewWords() {
  try {
    if (heldNewWordIds.size === 0) {
      localStorage.removeItem(
        HELD_NEW_WORDS_STORAGE_KEY
      );
      return;
    }

    localStorage.setItem(
      HELD_NEW_WORDS_STORAGE_KEY,
      JSON.stringify([...heldNewWordIds])
    );
  } catch (error) {
    console.error(
      "Erro ao salvar palavras guardadas:",
      error
    );
  }
}

function pruneHeldNewWords() {
  const validIds = new Set(
    allCards.map((card) => String(card.id))
  );

  const previousSize = heldNewWordIds.size;

  heldNewWordIds = new Set(
    [...heldNewWordIds].filter((id) =>
      validIds.has(id)
    )
  );

  if (heldNewWordIds.size !== previousSize) {
    saveHeldNewWords();
  }

  updateHeldNewWordsUI();
}

function updateHeldNewWordsUI() {
  restoreHeldWordsButton.textContent =
    `Recolocar (${heldNewWordIds.size})`;

  restoreHeldWordsButton.disabled =
    heldNewWordIds.size === 0 ||
    isManagingNewWordsDeck;

  const hasCurrentCard =
    Boolean(cards[currentIndex]);

  holdCurrentWordButton.disabled =
    !isNewWordsMode ||
    !hasCurrentCard ||
    isChangingCard ||
    writeComparisonDone ||
    isManagingNewWordsDeck;
}

async function loadCards() {
  try {
    const response = await fetch(DATA_URL);

    if (!response.ok) {
      throw new Error("Não foi possível carregar o arquivo words.json.");
    }

    const data = await response.json();

    allCards = data
      .filter((card) => card.active)
      .sort((a, b) => a.order - b.order);

    pruneHeldNewWords();

    buildExerciseVocabularyIndex();

    stats = loadStats();

    if (allCards.length === 0) {
      setupMessage.textContent = "Nenhum card ativo encontrado no JSON.";
      return;
    }

    fillFilterOptions();
    updateModeUI();
  } catch (error) {
    console.error(error);
    setupMessage.textContent = "Erro ao carregar os cards. Verifique data/words.json.";
  }
}

function loadStats() {
  const savedStats = localStorage.getItem(STORAGE_KEY);

  if (!savedStats) {
    return {};
  }

  try {
    return JSON.parse(savedStats);
  } catch (error) {
    console.error("Erro ao ler estatísticas salvas:", error);
    return {};
  }
}

function saveStats() {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(stats));
}

function getCardStats(cardId) {
  const key = String(cardId);

  if (!stats[key]) {
    stats[key] = {
      correct: 0,
      wrong: 0,
      seen: 0,
      lastResult: null,
      lastSeen: null,

      // SRS imediato para treino contínuo
      mastery: 0,
      streak: 0,
      lapses: 0,
      priority: 100,
      lastAnsweredAt: null
    };
  }

  // Migração para usuários que já têm estatísticas salvas antigas
  if (stats[key].mastery === undefined) {
    stats[key].mastery = calculateInitialMastery(stats[key]);
  }

  if (stats[key].streak === undefined) {
    stats[key].streak = 0;
  }

  if (stats[key].lapses === undefined) {
    stats[key].lapses = 0;
  }

  if (stats[key].priority === undefined) {
    stats[key].priority = calculateCardPriority(stats[key]);
  }

  if (stats[key].lastAnsweredAt === undefined) {
    stats[key].lastAnsweredAt = stats[key].lastSeen || null;
  }

  return stats[key];
}

function calculateInitialMastery(cardStats) {
  const correct = cardStats.correct || 0;
  const wrong = cardStats.wrong || 0;

  if (correct + wrong === 0) {
    return 0;
  }

  return clamp(correct - wrong * 2, -6, 8);
}

function updateCardStats(card, isCorrect) {
  const cardStats = getCardStats(card.id);
  const now = new Date();

  // Guarda o estado ANTES da resposta atual
  const previousMastery = cardStats.mastery || 0;
  const previousStreak = cardStats.streak || 0;

  cardStats.seen++;
  cardStats.lastSeen = now.toISOString();
  cardStats.lastAnsweredAt = now.toISOString();

  if (isCorrect) {
    cardStats.correct++;
    cardStats.lastResult = "correct";
    cardStats.streak = previousStreak + 1;

    cardStats.mastery = clamp(
      previousMastery + getCorrectMasteryGain(cardStats),
      -6,
      10
    );
  } else {
    cardStats.wrong++;
    cardStats.lastResult = "wrong";

    // Só é lapse se a palavra já estava aprendida
    const wasLearned =
      previousStreak >= 3 ||
      previousMastery >= 4;

    if (wasLearned) {
      cardStats.lapses = (cardStats.lapses || 0) + 1;
    }

    cardStats.streak = 0;

    cardStats.mastery = clamp(
      previousMastery - 3,
      -6,
      10
    );
  }

  cardStats.priority = calculateCardPriority(cardStats);

  saveStats();
}

function getCorrectMasteryGain(cardStats) {
  const streak = cardStats.streak || 0;

  if (streak >= 4) {
    return 0.75;
  }

  if (streak >= 2) {
    return 1;
  }

  return 1.5;
}

function getMaxSeenAmongCards(cardList = baseSessionCards) {
  if (!cardList || cardList.length === 0) {
    return 0;
  }

  return Math.max(
    0,
    ...cardList.map((card) => {
      return stats[String(card.id)]?.seen || 0;
    })
  );
}

function calculateCardPriority(cardStats) {
  const seen = cardStats.seen || 0;
  const wrong = cardStats.wrong || 0;
  const mastery = cardStats.mastery || 0;
  const lapses = cardStats.lapses || 0;

  let priority = 80;

  // Prioridade pela dificuldade
  priority += wrong * 14;
  priority += lapses * 6;
  priority -= mastery * 9;
  priority -= (cardStats.streak || 0) * 5;

  if (cardStats.lastResult === "wrong") {
    priority += 45;
  }

  // Prioridade por pouca exposição
  const maxSeen = getMaxSeenAmongCards();

  const exposureGap = clamp(
    maxSeen - seen,
    0,
    UNDEREXPOSURE_MAX_GAP
  );

  priority += exposureGap * UNDEREXPOSURE_PRIORITY_PER_VIEW;

  // Cards nunca vistos continuam com prioridade alta
  if (seen === 0) {
    priority += 50;
  }

  return clamp(Math.round(priority), 5, 220);
}

function fillFilterOptions() {
  fillCheckboxGroup(
    typeFilterGroup,
    getUniqueValues(allCards, (card) => card.grammar.type)
  );

  fillCheckboxGroup(
    themeFilterGroup,
    getUniqueThemeValues(allCards)
  );

  fillCheckboxGroup(
    sourceFilterGroup,
    getUniqueValues(allCards, (card) => getCardSource(card))
  );

  updateSourceSpecificFilters();
}

function fillCheckboxGroup(container, values) {
  container.innerHTML = "";

  values.forEach((value) => {
    const label = document.createElement("label");
    label.className = "checkbox-option";

    const checkbox = document.createElement("input");
    checkbox.type = "checkbox";
    checkbox.value = String(value);

    const text = document.createElement("span");
    text.textContent = String(value);

    label.append(checkbox, text);
    container.appendChild(label);
  });
}

function getCheckedValues(container) {
  return [...container.querySelectorAll('input[type="checkbox"]:checked')]
    .map((checkbox) => checkbox.value);
}

function hasCheckedValue(container, value) {
  return getCheckedValues(container).includes(value);
}

function getCardSource(card) {
  const classification = card.classification || {};

  if (classification.source) {
    return classification.source;
  }

  if (classification.sourceTitle || classification.sourceArtist) {
    return "música";
  }

  return "livro";
}

function getCardSourceTitle(card) {
  return card.classification?.sourceTitle || null;
}

function updateSourceSpecificFilters() {
  const selectedSources = getCheckedValues(sourceFilterGroup);

  const shouldShowChapters = selectedSources.includes("livro");
  const shouldShowSourceTitles = selectedSources.includes("música");

  chapterFilterLabel.classList.toggle("hidden", !shouldShowChapters);
  sourceTitleFilterLabel.classList.toggle("hidden", !shouldShowSourceTitles);

  if (!shouldShowChapters) {
    chapterFilterLabel.removeAttribute("open");
  }

  if (!shouldShowSourceTitles) {
    sourceTitleFilterLabel.removeAttribute("open");
  }

  if (shouldShowChapters) {
    fillCheckboxGroup(
      chapterFilterGroup,
      getUniqueValues(
        allCards.filter((card) => getCardSource(card) === "livro"),
        (card) => card.classification.chapter
      )
    );
  } else {
    chapterFilterGroup.innerHTML = "";
  }

  if (shouldShowSourceTitles) {
    fillCheckboxGroup(
      sourceTitleFilterGroup,
      getUniqueValues(
        allCards.filter((card) => getCardSource(card) === "música"),
        (card) => getCardSourceTitle(card)
      )
    );
  } else {
    sourceTitleFilterGroup.innerHTML = "";
  }
}

function getUniqueThemeValues(cardList) {
  return [...new Set(cardList.flatMap(getCardThemes))]
    .sort((a, b) => String(a).localeCompare(String(b), "pt-BR"));
}

function getCardThemes(card) {
  const classification = card.classification || {};

  if (Array.isArray(classification.themes)) {
    return classification.themes.filter(Boolean);
  }

  if (Array.isArray(classification.theme)) {
    return classification.theme.filter(Boolean);
  }

  if (classification.theme) {
    return [classification.theme];
  }

  return [];
}

function getUniqueValues(list, getter) {
  return [...new Set(list.map(getter).filter((value) => value !== null && value !== undefined))]
    .sort((a, b) => {
      if (typeof a === "number" && typeof b === "number") {
        return a - b;
      }

      return String(a).localeCompare(String(b), "pt-BR");
    });
}

function getFilteredCards() {
  const selectedTypes = getCheckedValues(typeFilterGroup);
  const selectedThemes = getCheckedValues(themeFilterGroup);
  const selectedSources = getCheckedValues(sourceFilterGroup);
  const selectedChapters = getCheckedValues(chapterFilterGroup);
  const selectedSourceTitles = getCheckedValues(sourceTitleFilterGroup);

  return allCards.filter((card) => {
    const matchesType =
      selectedTypes.length === 0 ||
      selectedTypes.includes(card.grammar.type);

    const cardThemes = getCardThemes(card);

    const matchesTheme =
      selectedThemes.length === 0 ||
      cardThemes.some((theme) => selectedThemes.includes(theme));

    const cardSource = getCardSource(card);

    const matchesSource =
      selectedSources.length === 0 ||
      selectedSources.includes(cardSource);

    const matchesChapter =
      cardSource !== "livro" ||
      selectedChapters.length === 0 ||
      selectedChapters.includes(String(card.classification.chapter));

    const matchesSourceTitle =
      cardSource !== "música" ||
      selectedSourceTitles.length === 0 ||
      selectedSourceTitles.includes(getCardSourceTitle(card));

    return (
      matchesType &&
      matchesTheme &&
      matchesSource &&
      matchesChapter &&
      matchesSourceTitle
    );
  });
}

function startSession() {
  baseSessionCards = getCardsAvailableForDirection(getFilteredCards());

  if (baseSessionCards.length === 0) {
    setupMessage.textContent = "Nenhum card encontrado com esses filtros e essa direção.";
    return;
  }

  if (isNewWordsMode) {
    startNewWordsSession(baseSessionCards);
  } else {
    cards = buildFirstRoundDeck(baseSessionCards);
  }

  setupMessage.textContent = "";

  currentIndex = 0;
  correctCount = 0;
  wrongCount = 0;
  cardsSinceUnseen = 0;
  sessionAnswers = [];
  answerVisible = false;
  isChangingCard = false;

  showScreen(studyScreen);

  startStudyTimer();
  showCard();
}

function getCardsAvailableForDirection(cardList) {
  const direction = directionSelect.value;

  if (direction === "audio-pt" || direction === "audio-sv") {
    return cardList.filter((card) => card.media?.audio?.src);
  }

  return cardList;
}

function startNewWordsSession(filteredCards) {
  const availableCards = filteredCards.filter(
    (card) =>
      !heldNewWordIds.has(String(card.id))
  );

  newWordsPool = getNewWordsPool(
    availableCards
  );

  newWordsActiveCount =
    newWordsPool.length > 0 ? 1 : 0;

  newWordsRoundResults = {};
  newWordsCardProgress = {};

  allNewWordsIntroduced =
    newWordsPool.length <= 1;

  cards = newWordsPool
    .slice(0, newWordsActiveCount)
    .map(createStudyOccurrence);
}

function getNewWordProgress(cardId) {
  const key = String(cardId);

  if (!newWordsCardProgress[key]) {
    newWordsCardProgress[key] = {
      attempts: 0,
      hadWrongAnswer: false,
      awaitingConfirmation: false,
      confirmed: false,
      delayRounds: 0
    };
  }

  return newWordsCardProgress[key];
}

function updateNewWordProgress(card, isCorrect) {
  const progress = getNewWordProgress(card.id);

  progress.attempts++;

  if (!isCorrect) {
    // Se errou alguma vez, perde o benefício de ter acertado de primeira.
    progress.hadWrongAnswer = true;
    progress.awaitingConfirmation = false;
    progress.confirmed = false;
    progress.delayRounds = 0;
    return;
  }

  if (
    progress.attempts === 1 &&
    !progress.hadWrongAnswer
  ) {
    // Acertou a palavra logo na primeira apresentação.
    progress.awaitingConfirmation = true;
    progress.confirmed = false;
    progress.delayRounds = NEW_WORD_FIRST_TRY_DELAY;
    return;
  }

  if (progress.awaitingConfirmation) {
    // Acertou novamente depois do período de espera.
    progress.awaitingConfirmation = false;
    progress.confirmed = true;
    progress.delayRounds = 0;
  }
}

function advanceNewWordsConfirmationDelays() {
  Object.values(newWordsCardProgress).forEach((progress) => {
    if (
      progress.awaitingConfirmation &&
      progress.delayRounds > 0
    ) {
      progress.delayRounds--;
    }
  });
}

function canNewWordAdvance(card) {
  const progress = getNewWordProgress(card.id);

  // Uma palavra acertada de primeira não impede a entrada
  // de novas palavras enquanto espera pela confirmação.
  if (progress.awaitingConfirmation || progress.confirmed) {
    return true;
  }

  return newWordsRoundResults[String(card.id)] === true;
}

function shouldShowNewWord(card) {
  const progress = getNewWordProgress(card.id);

  if (progress.confirmed) {
    return false;
  }

  if (
    progress.awaitingConfirmation &&
    progress.delayRounds > 0
  ) {
    return false;
  }

  return true;
}

function buildNewWordsRoundDeck(activeCards) {
  let nextCards = activeCards.filter(shouldShowNewWord);

  /*
   * Pode acontecer de todas as palavras ativas estarem esperando
   * confirmação. Nesse caso, mostramos a palavra que estiver mais
   * próxima de terminar seu período de espera, evitando um baralho vazio.
   */
  if (nextCards.length === 0) {
    const awaitingCards = activeCards
      .filter((card) => {
        const progress = getNewWordProgress(card.id);
        return progress.awaitingConfirmation;
      })
      .sort((cardA, cardB) => {
        const progressA = getNewWordProgress(cardA.id);
        const progressB = getNewWordProgress(cardB.id);

        return progressA.delayRounds - progressB.delayRounds;
      });

    if (awaitingCards.length > 0) {
      const confirmationCard = awaitingCards[0];
      const progress = getNewWordProgress(confirmationCard.id);

      progress.delayRounds = 0;
      nextCards = [confirmationCard];
    }
  }

  /*
   * Depois que todas tiverem sido confirmadas, mantém o comportamento
   * contínuo antigo com o SRS geral.
   */
  if (nextCards.length === 0) {
    return buildWeightedDeck(activeCards);
  }

  return shuffleArray(nextCards.map(createStudyOccurrence));
}

function getNewWordsPool(filteredCards) {
  return shuffleArray(filteredCards);
}

function clamp(value, min, max) {
  return Math.min(Math.max(value, min), max);
}

function shuffleArray(array) {
  const copy = [...array];

  for (let i = copy.length - 1; i > 0; i--) {
    const randomIndex = Math.floor(Math.random() * (i + 1));
    [copy[i], copy[randomIndex]] = [copy[randomIndex], copy[i]];
  }

  return copy;
}

function buildFirstRoundDeck(cardList) {
  return buildSrsDeck(cardList);
}

function openGrammarScreen() {
  stopStudyTimer();
  setupMessage.textContent = "";

  showScreen(grammarScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function backFromGrammarScreen() {
  showScreen(setupScreen);
}

function openPronunciationRules() {
  stopStudyTimer();
  setupMessage.textContent = "";

  showScreen(pronunciationScreen);

  renderPronunciationRules();

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function backFromPronunciationRules() {
  showScreen(grammarScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function openPluralRules() {
  stopStudyTimer();
  setupMessage.textContent = "";

  showScreen(pluralScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function backFromPluralRules() {
  showScreen(grammarScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function openVerbForms() {
  stopStudyTimer();
  setupMessage.textContent = "";

  showScreen(verbFormsScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function backFromVerbForms() {
  showScreen(grammarScreen);

  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function renderPronunciationRules() {
  if (!pronunciationRulesList) {
    return;
  }

  pronunciationRulesList.innerHTML = "";

  PRONUNCIATION_RULES.forEach((rule) => {
    const article = document.createElement("article");
    article.className = "pronunciation-rule";

    const title = document.createElement("h3");
    title.textContent = rule.title;

    const text = document.createElement("p");
    text.textContent = rule.text;

    const examples = document.createElement("div");
    examples.className = "pronunciation-examples";

    rule.examples
      .map(findCardWithAudioBySwedishTerm)
      .filter(Boolean)
      .forEach((card) => {
        examples.appendChild(createPronunciationExampleButton(card));
      });

    if (examples.children.length === 0) {
      const empty = document.createElement("p");
      empty.className = "pronunciation-empty";
      empty.textContent = "Nenhum exemplo com áudio encontrado no words.json para esta regra.";
      examples.appendChild(empty);
    }

    article.append(title, text, examples);
    pronunciationRulesList.appendChild(article);
  });
}

function findCardWithAudioBySwedishTerm(term) {
  const normalizedTerm = normalizePronunciationTerm(term);

  return allCards.find((card) => {
    const swedish = card.term?.swedish || "";
    const audioSrc = card.media?.audio?.src;

    return audioSrc && normalizePronunciationTerm(swedish) === normalizedTerm;
  });
}

function normalizePronunciationTerm(value) {
  return String(value)
    .trim()
    .toLocaleLowerCase("sv-SE")
    .replace(/[.!?]+$/g, "");
}

function createPronunciationExampleButton(card) {
  const button = document.createElement("button");
  button.type = "button";
  button.className = "pronunciation-example";
  button.dataset.audioSrc = card.media.audio.src;

  const swedish = document.createElement("strong");
  swedish.textContent = card.term.swedish;

  const portuguese = document.createElement("span");
  portuguese.textContent = card.term.portuguese;

  const action = document.createElement("em");
  action.textContent = "Ouvir";

  button.append(swedish, portuguese, action);
  return button;
}

function playPronunciationExample(event) {
  const button = event.target.closest(".pronunciation-example");

  if (!button) {
    return;
  }

  const audioSrc = button.dataset.audioSrc;

  if (!audioSrc) {
    return;
  }

  const audio = new Audio(audioSrc);
  audio.play().catch((error) => {
    console.error("Erro ao tocar áudio do exemplo:", error);
  });
}

function backToSetup() {
  stopStudyTimer();
  showScreen(setupScreen);

  flashcard.classList.remove("flipped", "correct-preview", "wrong-preview");
  answerVisible = false;
  isChangingCard = false;
}

function repeatSession() {
  if (baseSessionCards.length === 0) {
    backToSetup();
    return;
  }

  if (isNewWordsMode) {
    startNewWordsSession(baseSessionCards);
  } else {
    cards = buildFirstRoundDeck(baseSessionCards);
  }

  currentIndex = 0;
  correctCount = 0;
  wrongCount = 0;
  cardsSinceUnseen = 0;
  sessionAnswers = [];
  answerVisible = false;
  isChangingCard = false;

  showScreen(studyScreen);

  startStudyTimer();
  showCard();
}

function showEmptyDeckState() {
  currentIndex = 0;
  answerVisible = false;
  writeComparisonDone = false;

  flashcard.classList.remove(
    "flipped",
    "correct-preview",
    "wrong-preview"
  );

  hideQuestionAudio();

  pronunciationButton.classList.add("hidden");

  audioSwedishText.textContent = "";
  audioSwedishText.classList.add("hidden");

  showSwedishFromAudioButton.classList.add(
    "hidden"
  );

  showSwedishFromAudioButton.disabled = true;

  questionLabel.textContent = "Palavras novas";

  updateQuestionText(
    "Nenhuma palavra disponível"
  );

  answerLabel.textContent = "";
  answerText.textContent = "";

  answerGrammarForm.textContent = "";
  answerGrammarForm.classList.add("hidden");

  writeBox.classList.add("hidden");
  resultButtons.classList.add("hidden");

  correctButton.disabled = true;
  wrongButton.disabled = true;
  checkAnswerButton.disabled = true;
  nextWriteCardButton.disabled = true;

  backHint.classList.add("hidden");

  flashcard.setAttribute(
    "aria-disabled",
    "true"
  );

  if (heldNewWordIds.size > 0) {
    message.textContent =
      "Não há palavras disponíveis. Use Recolocar para devolver palavras ao deck.";
  } else {
    message.textContent =
      "Não há cards disponíveis nesta sessão.";
  }

  updateHeldNewWordsUI();
}

function showCard() {
  const card = cards[currentIndex];

  if (!card) {
    showEmptyDeckState();
    return;
  }

  message.textContent = "";

  const content = getCardContent(card);

  answerVisible = false;

  flashcard.classList.remove("flipped", "correct-preview", "wrong-preview");

  updateQuestionLabel(content);
  updateQuestionText(content.question);
  answerText.textContent = content.answer;
  answerLabel.textContent = content.answerLabel;

  updateAnswerGrammarForm(card);

  showQuestionMedia(card, content.questionType);
  resetAudioSwedishHint(card, content);
  resetWriteMode();
  updateModeUI();

  correctButton.disabled = true;
  wrongButton.disabled = true;

  preloadLikelyAudios();

  updateHeldNewWordsUI();

}

async function holdCurrentNewWord() {
  if (
    !isNewWordsMode ||
    isChangingCard ||
    writeComparisonDone ||
    isManagingNewWordsDeck
  ) {
    return;
  }

  const currentCard = cards[currentIndex];

  if (!currentCard) {
    return;
  }

  isManagingNewWordsDeck = true;
  updateHeldNewWordsUI();

  const cardId = String(currentCard.id);

  heldNewWordIds.add(cardId);
  saveHeldNewWords();

  /*
   * Remove a palavra do pool principal.
   */
  newWordsPool = newWordsPool.filter(
    (card) =>
      String(card.id) !== cardId
  );

  /*
   * Apaga somente o progresso temporário
   * do modo Palavras novas.
   * Não apaga as estatísticas gerais.
   */
  delete newWordsRoundResults[cardId];
  delete newWordsCardProgress[cardId];

  /*
   * Remove todas as ocorrências dessa
   * palavra da rodada, inclusive eventual
   * repetição agendada depois de um erro.
   */
  const cardsBeforeCurrent = cards
    .slice(0, currentIndex)
    .filter(
      (card) =>
        String(card.id) !== cardId
    );

  const cardsAfterCurrent = cards
    .slice(currentIndex + 1)
    .filter(
      (card) =>
        String(card.id) !== cardId
    );

  cards = [
    ...cardsBeforeCurrent,
    ...cardsAfterCurrent
  ];

  const nextIndex =
  cardsBeforeCurrent.length;

  const hasNextCardInCurrentRound =
    cardsAfterCurrent.length > 0;

  cardAudio.pause();

  /*
   * Se não sobrou nenhuma palavra.
   */
  if (newWordsPool.length === 0) {
    newWordsActiveCount = 0;
    allNewWordsIntroduced = true;

    cards = [];
    currentIndex = 0;

    showCard();
  } else {
    /*
     * Mantém a quantidade ativa dentro
     * dos limites do novo pool.
     */
    newWordsActiveCount = Math.min(
      Math.max(newWordsActiveCount, 1),
      newWordsPool.length
    );

    allNewWordsIntroduced =
      newWordsActiveCount >=
      newWordsPool.length;

    /*
     * Se havia outro card depois do atual,
     * ele ocupa automaticamente a posição
     * do card removido.
     */
    if (hasNextCardInCurrentRound) {
      currentIndex = nextIndex;
      showCard();
    } else {
      /*
       * O card removido era o último
       * disponível nesta rodada.
       * Inicia corretamente a próxima.
       */
      currentIndex = 0;

      updateNewWordsDeckAfterRound();

      showCard();
    }
  }

  if (cards.length === 0) {
    message.textContent =
      `Palavra guardada. ${heldNewWordIds.size} palavra(s) fora do deck. Use Recolocar para devolvê-la(s).`;
  } else {
    message.textContent =
      `Palavra guardada. ${heldNewWordIds.size} palavra(s) fora do deck.`;
  }

  await wait(180);

  isManagingNewWordsDeck = false;
  updateHeldNewWordsUI();
}

function restoreHeldNewWords() {
  if (
    heldNewWordIds.size === 0 ||
    isManagingNewWordsDeck
  ) {
    return;
  }

  isManagingNewWordsDeck = true;
  updateHeldNewWordsUI();

  const restoredIds =
    new Set(heldNewWordIds);

  /*
   * Libera todas as palavras guardadas.
   */
  heldNewWordIds.clear();
  saveHeldNewWords();

  /*
   * Evita duplicação no pool atual.
   */
  const existingPoolIds = new Set(
    newWordsPool.map(
      (card) => String(card.id)
    )
  );

  /*
   * Das palavras restauradas, só entram
   * imediatamente nesta sessão aquelas
   * que pertencem aos filtros atualmente
   * escolhidos.
   */
  const restoredCards =
    baseSessionCards.filter((card) => {
      const id = String(card.id);

      return (
        restoredIds.has(id) &&
        !existingPoolIds.has(id)
      );
    });

  /*
   * Acrescenta-as no fim do pool,
   * embaralhadas, sem destruir a ordem/
   * progresso das palavras que já estavam
   * sendo estudadas.
   */
  newWordsPool.push(
    ...shuffleArray(restoredCards)
  );

  if (
    newWordsPool.length > 0 &&
    newWordsActiveCount === 0
  ) {
    newWordsActiveCount = 1;
  }

  if (newWordsPool.length > 0) {
    newWordsActiveCount = Math.min(
      Math.max(newWordsActiveCount, 1),
      newWordsPool.length
    );
  }

  allNewWordsIntroduced =
    newWordsActiveCount >=
    newWordsPool.length;

  /*
   * Caso especial:
   * todas as cartas tinham sido guardadas.
   */
  if (
    cards.length === 0 &&
    newWordsPool.length > 0
  ) {
    currentIndex = 0;

    const activeCards =
      newWordsPool.slice(
        0,
        newWordsActiveCount
      );

    cards =
      buildNewWordsRoundDeck(activeCards);

    showCard();
  }

  isManagingNewWordsDeck = false;
  updateHeldNewWordsUI();

  if (restoredCards.length > 0) {
    message.textContent =
      `${restoredCards.length} palavra(s) recolocada(s) no deck.`;
  } else {
    message.textContent =
      "As palavras guardadas foram liberadas.";
  }
}

function preloadLikelyAudios() {
  const upcomingCards = cards.slice(currentIndex, currentIndex + 4);

  const priorityCards = [...baseSessionCards]
    .filter((card) => card?.media?.audio?.src)
    .sort((a, b) => {
      const aStats = getCardStats(a.id);
      const bStats = getCardStats(b.id);

      return calculateCardPriority(bStats) - calculateCardPriority(aStats);
    })
    .slice(0, 4);

  preloadAudioForCards([...upcomingCards, ...priorityCards], 8);
}

function resetAudioSwedishHint(card, content) {
  audioSwedishText.textContent = "";
  audioSwedishText.classList.add("hidden");

  showSwedishFromAudioButton.textContent = "Mostrar sueco";
  showSwedishFromAudioButton.classList.add("hidden");
  showSwedishFromAudioButton.disabled = true;

  if (content.questionType !== "audio") {
    return;
  }

  if (directionSelect.value === "audio-sv") {
    return;
  }

  showSwedishFromAudioButton.classList.remove("hidden");
  showSwedishFromAudioButton.disabled = false;
}

function showSwedishFromAudio() {
  const card = cards[currentIndex];

  if (!card) {
    return;
  }

  audioSwedishText.textContent = getStudySwedishText(card);
  audioSwedishText.classList.remove("hidden");

  showSwedishFromAudioButton.disabled = true;
  showSwedishFromAudioButton.textContent = "Sueco mostrado";
}

function updateQuestionText(text) {
  if (!text) {
    questionText.textContent = "";
    questionText.classList.add("hidden");
    return;
  }

  questionText.textContent = text;
  questionText.classList.remove("hidden");
}

function updateQuestionLabel(content) {
  if (content.questionType === "audio") {
    questionLabel.textContent = "";
    questionLabel.classList.add("hidden");
    return;
  }

  questionLabel.textContent = content.questionLabel;
  questionLabel.classList.remove("hidden");
}

function updateAnswerGrammarForm(card) {
  answerGrammarForm.textContent = "";
  answerGrammarForm.classList.add("hidden");
}

function isNounWithGender(card) {
  return (
    card.grammar?.type === "substantivo" &&
    (card.grammar?.gender === "en" || card.grammar?.gender === "ett")
  );
}

function getAvailableNounForms(card) {
  if (!isNounWithGender(card)) {
    return [];
  }

  const forms = ["singularIndefinite"];

  if (card.grammar?.definiteSingular) {
    forms.push("singularDefinite");
  }

  if (card.grammar?.plural) {
    forms.push("pluralIndefinite");
  }

  if (card.grammar?.definitePlural) {
    forms.push("pluralDefinite");
  }

  return forms;
}

function chooseRandomNounForm(card) {
  const forms = getAvailableNounForms(card);

  if (forms.length === 0) {
    return null;
  }

  return forms[Math.floor(Math.random() * forms.length)];
}

function isVerb(card) {
  return card.grammar?.type === "verbo";
}

function getAvailableVerbForms(card) {
  if (!isVerb(card)) {
    return [];
  }

  const forms = [];

  if (card.grammar?.infinitive && card.term?.portugueseInfinitive) {
    forms.push("infinitive");
  }

  if (card.grammar?.past && card.term?.portuguesePast) {
    forms.push("past");
  }

  if (card.grammar?.supine && card.term?.portugueseSupine) {
    forms.push("supine");
  }

  return forms;
}

function chooseRandomVerbForm(card) {
  const forms = getAvailableVerbForms(card);

  if (forms.length === 0) {
    return null;
  }

  return forms[Math.floor(Math.random() * forms.length)];
}

function createStudyOccurrence(card) {
  const occurrence = {
    ...card
  };

  const isTextDirection =
    directionSelect.value === "sv-pt" ||
    directionSelect.value === "pt-sv";

  // No modo "Novas palavras", sempre usar a forma básica
  if (isNewWordsMode) {
    if (isNounWithGender(card)) {
      occurrence._nounForm = "singularIndefinite";
    }

    if (isVerb(card)) {
      occurrence._verbForm = "infinitive";
    }

    return occurrence;
  }

  // Nos outros modos, mantém o comportamento atual
  if (isTextDirection && isNounWithGender(card)) {
    occurrence._nounForm = chooseRandomNounForm(card);
  }

  if (isTextDirection && isVerb(card)) {
    occurrence._verbForm = chooseRandomVerbForm(card);
  }

  return occurrence;
}

function getStudySwedishText(card) {
  const swedish = card.term?.swedish || "";

  if (isVerb(card) && card._verbForm) {
    switch (card._verbForm) {
      case "infinitive":
        return card.grammar.infinitive || swedish;

      case "present":
        return card.grammar.present || swedish;

      case "past":
        return card.grammar.past || swedish;

      case "supine":
        return card.grammar.supine || swedish;

      default:
        return swedish;
    }
  }

  if (!isNounWithGender(card)) {
    return swedish;
  }

  if (!card._nounForm) {
    return `${card.grammar.gender} ${swedish}`;
  }

  switch (card._nounForm) {
    case "singularIndefinite":
      return `${card.grammar.gender} ${swedish}`;

    case "singularDefinite":
      return card.grammar.definiteSingular || swedish;

    case "pluralIndefinite":
      return card.grammar.plural || swedish;

    case "pluralDefinite":
      return card.grammar.definitePlural || swedish;

    default:
      return `${card.grammar.gender} ${swedish}`;
  }
}

function getStudyPortugueseText(card) {
  const portuguese = card.term?.portuguese || "";

  if (isVerb(card) && card._verbForm) {
    switch (card._verbForm) {
      case "infinitive":
        return card.term.portugueseInfinitive || portuguese;

      case "present":
        return card.term.portuguesePresent || portuguese;

      case "past":
        return card.term.portuguesePast || portuguese;

      case "supine":
        return card.term.portugueseSupine || portuguese;

      default:
        return portuguese;
    }
  }

  if (!isNounWithGender(card) || !card._nounForm) {
    return portuguese;
  }

  const plural = card.term?.portuguesePlural || portuguese;
  const gender = card.term?.portugueseGender;

  const singularArticle =
    gender === "feminine" ? "a" : "o";

  const pluralArticle =
    gender === "feminine" ? "as" : "os";

  switch (card._nounForm) {
    case "singularIndefinite":
      return portuguese;

    case "singularDefinite":
      return `${singularArticle} ${portuguese}`;

    case "pluralIndefinite":
      return plural;

    case "pluralDefinite":
      return `${pluralArticle} ${plural}`;

    default:
      return portuguese;
  }
}

function getCardContent(card) {
  const direction = directionSelect.value;
  const studySwedish = getStudySwedishText(card);
  const studyPortuguese = getStudyPortugueseText(card);

  if (direction === "pt-sv") {
    return {
      questionType: "text",
      question: studyPortuguese,
      questionLabel: "Português",
      answer: studySwedish,
      answerLabel: "Sueco"
    };
  }

  if (direction === "audio-pt") {
    return {
      questionType: "audio",
      question: "",
      questionLabel: "Áudio em sueco",
      answer: card.term.portuguese,
      answerLabel: "Português"
    };
  }

  if (direction === "audio-sv") {
    return {
      questionType: "audio",
      question: "",
      questionLabel: "Áudio em sueco",
      answer: studySwedish,
      answerLabel: "Sueco"
    };
  }

  return {
    questionType: "text",
    question: studySwedish,
    questionLabel: "Sueco",
    answer: studyPortuguese,
    answerLabel: "Português"
  };
}

function showQuestionMedia(card, questionType) {
  hideQuestionAudio();

  updatePronunciationButton(card);

  if (questionType === "audio") {
    showAudioIfAvailable(card);
  }
}

function updatePronunciationButton(card) {
  const audioSrc = card?.media?.audio?.src;
  const direction = directionSelect.value;

  if (!audioSrc || (direction !== "sv-pt" && direction !== "pt-sv")) {
    pronunciationButton.classList.add("hidden");
    return;
  }

  const targetFace =
    direction === "sv-pt" ?
    flashcard.querySelector(".card-front") :
    flashcard.querySelector(".card-back");

  targetFace.appendChild(pronunciationButton);

  pronunciationButton.classList.remove("hidden");
}

function playCurrentCardPronunciation(event) {
  event.stopPropagation();

  const card = cards[currentIndex];
  const audioSrc = card?.media?.audio?.src;

  if (!audioSrc) {
    return;
  }

  cardAudio.pause();
  cardAudio.src = audioSrc;
  cardAudio.currentTime = 0;

  cardAudio.play().catch((error) => {
    console.error("Erro ao tocar pronúncia:", error);
  });
}

function showAudioIfAvailable(card) {
  const audio = card.media?.audio;

  if (!audio || !audio.src) {
    hideQuestionAudio();
    return;
  }

  cardAudio.src = audio.src;
  cardAudioWrap.classList.remove("hidden");
}

function hideQuestionAudio() {
  cardAudioWrap.classList.add("hidden");
  cardAudio.pause();
  cardAudio.removeAttribute("src");
  cardAudio.load();
}

function revealAnswer() {
  if (!cards[currentIndex]) {
    return;
  }

  if (isChangingCard) {
    return;
  }

  if (getAnswerMode() === "write") {
    return;
  }

  answerVisible = !answerVisible;

  flashcard.classList.toggle("flipped", answerVisible);

  correctButton.disabled = !answerVisible;
  wrongButton.disabled = !answerVisible;
}

function wait(ms) {
  return new Promise((resolve) => {
    setTimeout(resolve, ms);
  });
}

function registerCurrentAnswer(isCorrect) {
  const currentCard = cards[currentIndex];
  const wasUnseen = (getCardStats(currentCard.id).seen || 0) === 0;

  if (isCorrect) {
    correctCount++;
  } else {
    wrongCount++;
  }

  if (isNewWordsMode) {
    newWordsRoundResults[String(currentCard.id)] = isCorrect;
    updateNewWordProgress(currentCard, isCorrect);
  }

  sessionAnswers.push({
    cardId: currentCard.id,
    nounForm: currentCard._nounForm || null,
    verbForm: currentCard._verbForm || null,
    swedish: getStudySwedishText(currentCard),
    portuguese: getStudyPortugueseText(currentCard),
    isCorrect
  });

  updateCardStats(currentCard, isCorrect);

  if (!isNewWordsMode) {
    if (wasUnseen) {
      cardsSinceUnseen = 0;
    } else {
      cardsSinceUnseen++;
    }
  }

  if (!isCorrect) {
    scheduleImmediateRetry(currentCard);
    preloadAudioForCard(currentCard);
  }

  return currentCard;
}

function scheduleImmediateRetry(card) {
  const insertIndex = Math.min(
    cards.length,
    currentIndex + getImmediateRetryGap()
  );

  const alreadyScheduledSoon = cards
    .slice(currentIndex + 1, insertIndex + 1)
    .some((queuedCard) => isSameStudyOccurrence(queuedCard, card));

  if (alreadyScheduledSoon) {
    return;
  }

  cards.splice(insertIndex, 0, card);
}

function isSameStudyOccurrence(cardA, cardB) {
  return (
    String(cardA.id) === String(cardB.id) &&
    (cardA._nounForm || null) === (cardB._nounForm || null) &&
    (cardA._verbForm || null) === (cardB._verbForm || null)
  );
}

function getImmediateRetryGap() {
  const remainingCards = cards.length - currentIndex - 1;

  if (remainingCards <= 2) {
    return 2;
  }

  return 3 + Math.floor(Math.random() * 3);
}

async function markAnswer(isCorrect) {
  if (!cards[currentIndex]) {
    return;
  }

  if (!answerVisible || isChangingCard) {
    return;
  }

  isChangingCard = true;

  correctButton.disabled = true;
  wrongButton.disabled = true;

  registerCurrentAnswer(isCorrect);

  flashcard.classList.remove("flipped", "correct-preview", "wrong-preview");

  await wait(220);

  goToNextCard();

  isChangingCard = false;

  updateHeldNewWordsUI();
}

function getRecentAccuracy() {
  const recentAnswers = sessionAnswers.slice(
    -UNSEEN_ACCURACY_WINDOW
  );

  // Ainda não há respostas suficientes
  // para avaliar o desempenho recente.
  if (recentAnswers.length < UNSEEN_MIN_ANSWERS) {
    return null;
  }

  const correctAnswers = recentAnswers.filter((answer) => {
    return answer.isCorrect;
  }).length;

  return correctAnswers / recentAnswers.length;
}

function canIntroduceUnseenCard() {
  if (isNewWordsMode) {
    return false;
  }

  const accuracy = getRecentAccuracy();

  if (accuracy === null) {
    return false;
  }

  // Regra normal:
  // após 4 cards, libera uma palavra nova
  // se a taxa recente for de pelo menos 60%.
  if (
    cardsSinceUnseen >= UNSEEN_CARD_MAX_GAP &&
    accuracy >= UNSEEN_MIN_ACCURACY
  ) {
    return true;
  }

  // Regra de espera prolongada:
  // após 15 cards, reduz a exigência para 50%.
  if (
    cardsSinceUnseen >= UNSEEN_HARD_MAX_GAP &&
    accuracy >= UNSEEN_HARD_MIN_ACCURACY
  ) {
    return true;
  }

  return false;
}

function prioritizeUnseenCardIfNeeded() {
  if (!canIntroduceUnseenCard()) {
    return;
  }

  const unseenIndex = cards.findIndex((card, index) => {
    if (index < currentIndex) {
      return false;
    }

    const cardStats = getCardStats(card.id);

    return (cardStats.seen || 0) === 0;
  });

  // Não existem mais palavras nunca vistas nessa fila.
  if (unseenIndex === -1) {
    return;
  }

  // A próxima já é uma palavra nunca vista.
  if (unseenIndex === currentIndex) {
    return;
  }

  // Remove da posição atual...
  const [unseenCard] = cards.splice(unseenIndex, 1);

  // ...e coloca como próximo card.
  cards.splice(currentIndex, 0, unseenCard);
}

function goToNextCard() {
  if (cards.length === 0) {
    currentIndex = 0;
    showCard();
    return;
  }

  currentIndex++;

  if (currentIndex >= cards.length) {
    currentIndex = 0;

    if (isNewWordsMode) {
      updateNewWordsDeckAfterRound();
    } else {
      cards = buildWeightedDeck(baseSessionCards);
    }
  }

  if (!isNewWordsMode) {
    prioritizeUnseenCardIfNeeded();
  }

  showCard();
}

function updateNewWordsDeckAfterRound() {
  let activeCards = newWordsPool.slice(0, newWordsActiveCount);

  const allActiveCardsCanAdvance = activeCards.every((card) => {
    return canNewWordAdvance(card);
  });

  if (
    allActiveCardsCanAdvance &&
    newWordsActiveCount < newWordsPool.length
  ) {
    newWordsActiveCount++;
  }

  allNewWordsIntroduced =
    newWordsActiveCount >= newWordsPool.length;

  // Reduz uma unidade do tempo de espera das palavras adiadas.
  advanceNewWordsConfirmationDelays();

  newWordsRoundResults = {};

  activeCards = newWordsPool.slice(0, newWordsActiveCount);
  cards = buildNewWordsRoundDeck(activeCards);
}

function buildWeightedDeck(filteredCards) {
  return buildSrsDeck(filteredCards);
}

function buildSrsDeck(filteredCards) {
  if (filteredCards.length === 0) {
    return [];
  }

  const sortedCards = [...filteredCards].sort((cardA, cardB) => {
    const statsA = getCardStats(cardA.id);
    const statsB = getCardStats(cardB.id);

    const priorityA = calculateCardPriority(statsA);
    const priorityB = calculateCardPriority(statsB);

    // Primeiro: maior prioridade
    if (priorityA !== priorityB) {
      return priorityB - priorityA;
    }

    // Em empate: menos vista primeiro
    return (statsA.seen || 0) - (statsB.seen || 0);
  });

  const firstPass = sortedCards.map(createStudyOccurrence);
  const extraCards = [];

  sortedCards.forEach((card) => {
    const cardStats = getCardStats(card.id);
    const repeatCount = getImmediateSrsRepeatCount(cardStats);

    // A primeira ocorrência já está em firstPass.
    for (let i = 1; i < repeatCount; i++) {
      extraCards.push(createStudyOccurrence(card));
    }
  });

  return [
    ...shuffleArray(firstPass),
    ...shuffleArray(extraCards)
  ];
}

function getImmediateSrsRepeatCount(cardStats) {
  const priority = calculateCardPriority(cardStats);

  if (priority >= 150) {
    return 5;
  }

  if (priority >= 115) {
    return 4;
  }

  if (priority >= 80) {
    return 3;
  }

  if (priority >= 45) {
    return 2;
  }

  return 1;
}

function showSummary() {
  stopStudyTimer();
  showScreen(summaryScreen);

  const total = correctCount + wrongCount;
  const correctRate = total === 0 ? 0 : Math.round((correctCount / total) * 100);

  summaryTotal.textContent = total;
  summaryCorrect.textContent = correctCount;
  summaryWrong.textContent = wrongCount;
  summaryCorrectRate.textContent = `${correctRate}%`;
  summaryTime.textContent = formatTime(getElapsedStudySeconds());

  renderWrongList();
}

function renderWrongList() {
  summaryWrongList.innerHTML = "";

  const wrongMap = new Map();

  sessionAnswers
    .filter((answer) => !answer.isCorrect)
    .forEach((answer) => {
      const occurrenceKey = [
        answer.cardId,
        answer.nounForm || "",
        answer.verbForm || ""
      ].join("|");

      const current = wrongMap.get(occurrenceKey) || {
        swedish: answer.swedish,
        portuguese: answer.portuguese,
        count: 0
      };

      current.count++;
      wrongMap.set(occurrenceKey, current);
    });

  const wrongItems = [...wrongMap.values()]
    .sort((a, b) => b.count - a.count)
    .slice(0, 5);

  if (wrongItems.length === 0) {
    const item = document.createElement("li");
    item.textContent = "Nenhum erro nesta sessão.";
    summaryWrongList.appendChild(item);
    return;
  }

  wrongItems.forEach((item) => {
    const listItem = document.createElement("li");
    listItem.textContent = `${item.swedish} — ${item.portuguese} (${item.count} erro${item.count > 1 ? "s" : ""})`;
    summaryWrongList.appendChild(listItem);
  });
}

function showError(errorMessage) {
  questionText.textContent = "Erro";
  answerText.textContent = "";
  flashcard.classList.remove("flipped");
  message.textContent = errorMessage;

  correctButton.disabled = true;
  wrongButton.disabled = true;
}

function getAnswerMode() {
  return answerModeSelect.value;
}

function updateModeUI() {
  const canWriteInSwedish = isSwedishAnswerModeAvailable();
  const forceWriteMode = directionSelect.value === "audio-sv";

  if (forceWriteMode) {
    answerModeSelect.value = "write";
    answerModeLabel.classList.add("hidden");
  } else if (!canWriteInSwedish) {
    answerModeSelect.value = "think";
    answerModeLabel.classList.add("hidden");
  } else {
    answerModeLabel.classList.remove("hidden");
  }

  const shouldShowWriteBox =
    canWriteInSwedish && getAnswerMode() === "write";

  if (shouldShowWriteBox) {
    writeBox.classList.remove("hidden");
    resultButtons.classList.add("hidden");
    flashcard.setAttribute("aria-disabled", "true");
  } else {
    writeBox.classList.add("hidden");
    resultButtons.classList.remove("hidden");
    flashcard.removeAttribute("aria-disabled");
  }
}

function resetWriteMode() {
  writeComparisonDone = false;

  answerInput.value = "";
  answerInput.disabled = false;

  checkAnswerButton.disabled = false;

  nextWriteCardButton.classList.add("hidden");
  nextWriteCardButton.disabled = true;

  writeResultBox.classList.add("hidden");
  userAnswerText.textContent = "";

  comparisonFeedbackText.textContent = "";
  comparisonFeedbackText.classList.remove("correct", "wrong");

  backHint.classList.remove("hidden");
}

function checkWrittenAnswer() {
  if (isChangingCard || writeComparisonDone) {
    return;
  }

  if (!isSwedishAnswerModeAvailable()) {
    message.textContent = "Modo escrita só está disponível quando a resposta é em sueco.";
    return;
  }

  const card = cards[currentIndex];

  if (!card) {
    showEmptyDeckState();
    return;
  }

  const content = getCardContent(card);

  const userAnswer = answerInput.value;
  const expectedAnswer = content.answer;

  if (!userAnswer.trim()) {
    message.textContent = "Digite uma resposta antes de comparar.";
    answerInput.focus();
    return;
  }

  const comparison = compareWrittenAnswer(userAnswer, expectedAnswer, card);

  answerVisible = true;
  writeComparisonDone = true;

  flashcard.classList.add("flipped");
  flashcard.classList.toggle("correct-preview", comparison.isCorrect);
  flashcard.classList.toggle("wrong-preview", !comparison.isCorrect);

  userAnswerText.textContent = userAnswer;
  comparisonFeedbackText.textContent = comparison.feedback;
  comparisonFeedbackText.classList.toggle("correct", comparison.isCorrect);
  comparisonFeedbackText.classList.toggle("wrong", !comparison.isCorrect);
  writeResultBox.classList.remove("hidden");

  backHint.classList.add("hidden");

  answerInput.disabled = true;
  checkAnswerButton.disabled = true;

  nextWriteCardButton.classList.remove("hidden");
  nextWriteCardButton.disabled = false;

  registerCurrentAnswer(comparison.isCorrect);

  updateHeldNewWordsUI();

}

async function nextWrittenCard() {
  if (!writeComparisonDone || isChangingCard) {
    return;
  }

  isChangingCard = true;

  nextWriteCardButton.disabled = true;

  flashcard.classList.remove("flipped", "correct-preview", "wrong-preview");

  await wait(220);

  goToNextCard();

  isChangingCard = false;

  updateHeldNewWordsUI();
}

function compareWrittenAnswer(userAnswer, expectedAnswer, card) {
  const normalizedUserAnswer = normalizeAnswer(userAnswer);
  const normalizedExpectedAnswer = normalizeAnswer(expectedAnswer);

  const comparison = getDamerauLevenshteinComparison(
    normalizedUserAnswer,
    normalizedExpectedAnswer
  );

  const isCorrect = comparison.distance === 0;

  if (isCorrect) {
    return {
      isCorrect: true,
      distance: 0,
      feedback: "Resposta correta."
    };
  }

  const nounComparison = compareNounWithGenderAnswer(
    card,
    normalizedUserAnswer,
    normalizedExpectedAnswer
  );

  if (nounComparison) {
    return nounComparison;
  }

  return {
    isCorrect: false,
    distance: comparison.distance,
    feedback: buildComparisonFeedbackFromOperations(
      comparison.operations,
      normalizedUserAnswer,
      normalizedExpectedAnswer
    )
  };
}

function getDamerauLevenshteinComparison(source, target) {
  const sourceLength = source.length;
  const targetLength = target.length;

  const distances = Array.from({
      length: sourceLength + 1
    }, () =>
    Array(targetLength + 1).fill(0)
  );

  const steps = Array.from({
      length: sourceLength + 1
    }, () =>
    Array(targetLength + 1).fill(null)
  );

  for (let i = 0; i <= sourceLength; i++) {
    distances[i][0] = i;

    if (i > 0) {
      steps[i][0] = {
        type: "delete",
        sourceIndex: i - 1,
        targetIndex: 0
      };
    }
  }

  for (let j = 0; j <= targetLength; j++) {
    distances[0][j] = j;

    if (j > 0) {
      steps[0][j] = {
        type: "insert",
        sourceIndex: 0,
        targetIndex: j - 1
      };
    }
  }

  for (let i = 1; i <= sourceLength; i++) {
    for (let j = 1; j <= targetLength; j++) {
      const substitutionCost = source[i - 1] === target[j - 1] ? 0 : 1;

      let bestDistance = distances[i - 1][j - 1] + substitutionCost;
      let bestStep = {
        type: substitutionCost === 0 ? "match" : "substitute",
        sourceIndex: i - 1,
        targetIndex: j - 1
      };

      const deleteDistance = distances[i - 1][j] + 1;
      if (deleteDistance < bestDistance) {
        bestDistance = deleteDistance;
        bestStep = {
          type: "delete",
          sourceIndex: i - 1,
          targetIndex: j
        };
      }

      const insertDistance = distances[i][j - 1] + 1;
      if (insertDistance < bestDistance) {
        bestDistance = insertDistance;
        bestStep = {
          type: "insert",
          sourceIndex: i,
          targetIndex: j - 1
        };
      }

      if (
        i > 1 &&
        j > 1 &&
        source[i - 1] === target[j - 2] &&
        source[i - 2] === target[j - 1]
      ) {
        const transposeDistance = distances[i - 2][j - 2] + 1;

        if (transposeDistance < bestDistance) {
          bestDistance = transposeDistance;
          bestStep = {
            type: "transpose",
            sourceIndex: i - 2,
            targetIndex: j - 2
          };
        }
      }

      distances[i][j] = bestDistance;
      steps[i][j] = bestStep;
    }
  }

  const operations = [];
  let i = sourceLength;
  let j = targetLength;

  while (i > 0 || j > 0) {
    const step = steps[i][j];

    if (!step) {
      break;
    }

    if (step.type === "match") {
      i--;
      j--;
      continue;
    }

    if (step.type === "substitute") {
      operations.push({
        type: "substitute",
        userChar: source[i - 1],
        expectedChar: target[j - 1],
        position: j
      });

      i--;
      j--;
      continue;
    }

    if (step.type === "delete") {
      operations.push({
        type: "delete",
        userChar: source[i - 1],
        position: i
      });

      i--;
      continue;
    }

    if (step.type === "insert") {
      operations.push({
        type: "insert",
        expectedChar: target[j - 1],
        position: j
      });

      j--;
      continue;
    }

    if (step.type === "transpose") {
      operations.push({
        type: "transpose",
        firstChar: source[i - 2],
        secondChar: source[i - 1],
        position: j - 1
      });

      i -= 2;
      j -= 2;
      continue;
    }
  }

  operations.reverse();

  return {
    distance: distances[sourceLength][targetLength],
    operations
  };
}

function compareNounWithGenderAnswer(card, normalizedUserAnswer, normalizedExpectedAnswer) {
  if (!isNounWithGender(card)) {
    return null;
  }

  const expectedArticle = card.grammar.gender;
  const expectedWord = normalizeAnswer(card.term?.swedish || "");

  if (!expectedArticle || !expectedWord) {
    return null;
  }

  const userNounParts = splitNounAnswer(normalizedUserAnswer);

  if (!userNounParts.word) {
    return null;
  }

  const feedbackParts = [];
  let totalDistance = 0;

  if (!userNounParts.article) {
    feedbackParts.push("Faltou o artigo.");
    totalDistance += 1;
  } else if (userNounParts.article !== expectedArticle) {
    const articleComparison = getDamerauLevenshteinComparison(
      userNounParts.article,
      expectedArticle
    );

    totalDistance += articleComparison.distance;

    if (userNounParts.article === "en" || userNounParts.article === "ett") {
      feedbackParts.push(`O artigo correto é "${expectedArticle}".`);
    } else {
      feedbackParts.push(
        buildArticleFeedbackFromOperations(
          articleComparison.operations,
          userNounParts.article,
          expectedArticle
        )
      );
    }
  }

  if (userNounParts.word !== expectedWord) {
    const wordComparison = getDamerauLevenshteinComparison(
      userNounParts.word,
      expectedWord
    );

    totalDistance += wordComparison.distance;

    feedbackParts.push(
      buildNounWordFeedbackFromOperations(
        wordComparison.operations,
        userNounParts.word,
        expectedWord
      )
    );
  }

  const cleanFeedbackParts = feedbackParts.filter(Boolean);

  if (cleanFeedbackParts.length === 0) {
    return null;
  }

  if (totalDistance > 4) {
    return {
      isCorrect: false,
      distance: totalDistance,
      feedback: `A resposta ficou muito diferente. Leia a forma correta com atenção e tente memorizá-la: "${normalizedExpectedAnswer}".`
    };
  }

  return {
    isCorrect: false,
    distance: totalDistance,
    feedback: joinFeedbackParts(cleanFeedbackParts)
  };
}

function splitNounAnswer(normalizedAnswer) {
  const words = normalizedAnswer.split(/\s+/).filter(Boolean);

  if (words.length <= 1) {
    return {
      article: "",
      word: normalizedAnswer,
      hasArticleAttempt: false
    };
  }

  return {
    article: words[0],
    word: words.slice(1).join(" "),
    hasArticleAttempt: true
  };
}

function buildArticleFeedbackFromOperations(operations, userArticle, expectedArticle) {
  const relevantOperations = operations.filter((operation) => {
    return operation.type !== "match";
  });

  if (relevantOperations.length === 0) {
    return "";
  }

  return relevantOperations
    .map((operation) => {
      return formatArticleOperation(operation, userArticle, expectedArticle);
    })
    .filter(Boolean)
    .join(" ");
}

function formatArticleOperation(operation, userArticle, expectedArticle) {
  if (operation.type === "substitute") {
    const location = describeArticleCharacterLocation(expectedArticle, operation.position);

    return `${capitalizeFirstLetter(location)}, você escreveu "${operation.userChar}", mas o correto é "${operation.expectedChar}".`;
  }

  if (operation.type === "insert") {
    const location = describeArticleCharacterLocation(expectedArticle, operation.position);

    return `Faltou um "${operation.expectedChar}" ${location}.`;
  }

  if (operation.type === "delete") {
    const location = describeArticleCharacterLocation(userArticle, operation.position);

    return `Há um "${operation.userChar}" extra ${location}.`;
  }

  if (operation.type === "transpose") {
    const location = describeArticleCharacterLocation(expectedArticle, operation.position);

    return `Algumas letras parecem estar invertidas ${location}.`;
  }

  return "";
}

function describeArticleCharacterLocation(article, charPosition) {
  const index = clamp(charPosition - 1, 0, Math.max(article.length - 1, 0));
  const articlePart = getWordPart(index, article.length);

  return `${articlePart} do artigo`;
}

function joinFeedbackParts(parts) {
  if (parts.length === 1) {
    return parts[0];
  }

  if (parts.length === 2) {
    return `${removeFinalPeriod(parts[0])} e ${lowercaseFirstLetter(parts[1])}`;
  }

  const firstParts = parts.slice(0, -1).map(removeFinalPeriod);
  const lastPart = lowercaseFirstLetter(parts[parts.length - 1]);

  return `${firstParts.join(", ")} e ${lastPart}`;
}

function removeFinalPeriod(text) {
  return text.replace(/\.$/, "");
}

function lowercaseFirstLetter(text) {
  if (!text) {
    return text;
  }

  return text.charAt(0).toLowerCase() + text.slice(1);
}

function buildNounWordFeedbackFromOperations(operations, userWord, expectedWord) {
  const relevantOperations = operations.filter((operation) => {
    return operation.type !== "match";
  });

  if (relevantOperations.length === 0) {
    return "Resposta correta.";
  }

  return relevantOperations
    .map((operation) => {
      return formatNounWordOperation(operation, userWord, expectedWord);
    })
    .filter(Boolean)
    .join(" ");
}

function formatNounWordOperation(operation, userWord, expectedWord) {
  if (operation.type === "substitute") {
    const location = describeNounWordCharacterLocation(expectedWord, operation.position);

    return `${capitalizeFirstLetter(location)}, você escreveu "${operation.userChar}", mas o correto é "${operation.expectedChar}".`;
  }

  if (operation.type === "insert") {
    const location = describeNounWordCharacterLocation(expectedWord, operation.position);

    return `Faltou a letra "${operation.expectedChar}" ${location}.`;
  }

  if (operation.type === "delete") {
    const location = describeNounWordCharacterLocation(userWord, operation.position);

    return `Há uma letra extra "${operation.userChar}" ${location}.`;
  }

  if (operation.type === "transpose") {
    const location = describeNounWordCharacterLocation(expectedWord, operation.position);

    return `Algumas letras parecem estar invertidas ${location}.`;
  }

  return "";
}

function describeNounWordCharacterLocation(word, charPosition) {
  const index = clamp(charPosition - 1, 0, Math.max(word.length - 1, 0));
  const wordPart = getWordPart(index, word.length);

  return `${wordPart} da palavra`;
}

function buildComparisonFeedbackFromOperations(operations, userAnswer, expectedAnswer) {
  const relevantOperations = operations.filter((operation) => {
    return operation.type !== "match";
  });

  if (relevantOperations.length === 0) {
    return "Resposta correta.";
  }

  if (relevantOperations.length > 4) {
    return `A resposta ficou muito diferente. Leia a forma correta com atenção e tente memorizá-la: "${expectedAnswer}".`;
  }

  return relevantOperations
    .map((operation) => {
      return formatComparisonOperation(operation, userAnswer, expectedAnswer);
    })
    .filter(Boolean)
    .join(" ");
}

function formatComparisonOperation(operation, userAnswer, expectedAnswer) {
  if (operation.type === "substitute") {
    const location = describeCharacterLocation(expectedAnswer, operation.position);

    return `${capitalizeFirstLetter(location)}, você escreveu "${operation.userChar}", mas o correto é "${operation.expectedChar}".`;
  }

  if (operation.type === "insert") {
    const location = describeCharacterLocation(expectedAnswer, operation.position);

    return `Faltou a letra "${operation.expectedChar}" ${location}.`;
  }

  if (operation.type === "delete") {
    const location = describeCharacterLocation(userAnswer, operation.position);

    return `Há uma letra extra "${operation.userChar}" ${location}.`;
  }

  if (operation.type === "transpose") {
    const location = describeCharacterLocation(expectedAnswer, operation.position);

    return `Algumas letras parecem estar invertidas ${location}.`;
  }

  return "";
}

function describeCharacterLocation(text, charPosition) {
  const index = clamp(charPosition - 1, 0, Math.max(text.length - 1, 0));

  if (text[index] === " ") {
    return describeSpaceLocation(text, index);
  }

  const wordRanges = getWordRanges(text);
  const wordIndex = wordRanges.findIndex((range) => {
    return index >= range.start && index <= range.end;
  });

  if (wordIndex === -1) {
    return "em uma parte da resposta";
  }

  const wordRange = wordRanges[wordIndex];
  const wordLength = wordRange.end - wordRange.start + 1;
  const relativeIndex = index - wordRange.start;

  const wordPart = getWordPart(relativeIndex, wordLength);
  const ordinalWord = getOrdinalWord(wordIndex + 1);

  return `${wordPart} da ${ordinalWord} palavra`;
}

function describeSpaceLocation(text, spaceIndex) {
  const textBeforeSpace = text.slice(0, spaceIndex).trim();
  const textAfterSpace = text.slice(spaceIndex + 1).trim();

  const wordsBefore = textBeforeSpace ? textBeforeSpace.split(/\s+/).length : 0;
  const wordsAfter = textAfterSpace ? textAfterSpace.split(/\s+/).length : 0;

  if (wordsBefore > 0 && wordsAfter > 0) {
    const previousWord = getOrdinalWord(wordsBefore);
    const nextWord = getOrdinalWord(wordsBefore + 1);

    return `entre a ${previousWord} palavra e a ${nextWord} palavra`;
  }

  if (wordsBefore === 0) {
    return "antes da primeira palavra";
  }

  return `depois da ${getOrdinalWord(wordsBefore)} palavra`;
}

function getWordRanges(text) {
  const ranges = [];
  const wordRegex = /\S+/g;
  let match;

  while ((match = wordRegex.exec(text)) !== null) {
    ranges.push({
      start: match.index,
      end: match.index + match[0].length - 1,
      word: match[0]
    });
  }

  return ranges;
}

function getWordPart(relativeIndex, wordLength) {
  if (wordLength <= 1) {
    return "na palavra";
  }

  if (wordLength === 2) {
    return relativeIndex === 0 ? "no início" : "no final";
  }

  const firstThirdLimit = Math.ceil(wordLength / 3);
  const lastThirdStart = Math.floor((wordLength * 2) / 3);

  if (relativeIndex < firstThirdLimit) {
    return "no início";
  }

  if (relativeIndex >= lastThirdStart) {
    return "no final";
  }

  return "no meio";
}

function getOrdinalWord(number) {
  const ordinals = [
    "primeira",
    "segunda",
    "terceira",
    "quarta",
    "quinta",
    "sexta",
    "sétima",
    "oitava",
    "nona",
    "décima"
  ];

  return ordinals[number - 1] || `${number}ª`;
}

function isSwedishAnswerModeAvailable() {
  return (
    directionSelect.value === "pt-sv" ||
    directionSelect.value === "audio-sv"
  );
}

function normalizeAnswer(value) {
  return value
    .toLowerCase()
    .trim()
    .replace(/[.,!?;:()"]/g, "")
    .replace(/\s+/g, " ");
}

function applySavedTheme() {
  const savedTheme = localStorage.getItem(THEME_STORAGE_KEY);

  if (savedTheme === "dark") {
    document.body.classList.add("dark-theme");
  } else {
    document.body.classList.remove("dark-theme");
  }

  updateThemeButtonText();
}

function toggleTheme() {
  const isDark = document.body.classList.toggle("dark-theme");

  localStorage.setItem(THEME_STORAGE_KEY, isDark ? "dark" : "light");
  updateThemeButtonText();
}

function updateThemeButtonText() {
  const themeColorMeta = document.querySelector('meta[name="theme-color"]');

  if (document.body.classList.contains("dark-theme")) {
    themeToggleButton.textContent = "Tema claro";

    if (themeColorMeta) {
      themeColorMeta.setAttribute("content", "#101A2B");
    }
  } else {
    themeToggleButton.textContent = "Tema escuro";

    if (themeColorMeta) {
      themeColorMeta.setAttribute("content", "#FFF7D8");
    }
  }
}

function startStudyTimer() {
  sessionStartTime = Date.now();
  elapsedBeforePause = 0;

  if (timerIntervalId) {
    clearInterval(timerIntervalId);
    timerIntervalId = null;
  }
}

function stopStudyTimer() {
  if (timerIntervalId) {
    clearInterval(timerIntervalId);
    timerIntervalId = null;
  }
}

function getElapsedStudySeconds() {
  if (!sessionStartTime) {
    return 0;
  }

  return Math.floor((Date.now() - sessionStartTime) / 1000) + elapsedBeforePause;
}

function formatTime(totalSeconds) {
  const minutes = Math.floor(totalSeconds / 60);
  const seconds = totalSeconds % 60;

  return `${String(minutes).padStart(2, "0")}:${String(seconds).padStart(2, "0")}`;
}

function toggleNewWordsMode() {
  isNewWordsMode = !isNewWordsMode;
  updateNewWordsModeUI();
}

function updateNewWordsModeUI() {
  if (isNewWordsMode) {
    newWordsToggleButton.classList.add("active");
    newWordsToggleButton.textContent = "Palavras novas: ligado";
  } else {
    newWordsToggleButton.classList.remove("active");
    newWordsToggleButton.textContent = "Palavras novas";
  }
}

function capitalizeFirstLetter(text) {
  if (!text) {
    return "";
  }

  return text.charAt(0).toUpperCase() + text.slice(1);
}

function openWordsScreen(viewMode = "all") {
  stopStudyTimer();
  setupMessage.textContent = "";
  wordsViewMode = viewMode;

  expandedWordCardId = null;
  expandedWordsLetter = null;

  showScreen(wordsScreen);

  renderWordsList();
  window.scrollTo({
    top: 0,
    behavior: "smooth"
  });
}

function backFromWordsScreen() {
  showScreen(setupScreen);
}

function toggleWordsDirection() {
  wordsDirection = wordsDirection === "sv-pt" ? "pt-sv" : "sv-pt";
  expandedWordCardId = null;
  expandedWordsLetter = null;
  renderWordsList();
}

function renderWordsList() {
  if (!wordsList) {
    return;
  }

  const filteredCards = getCardsForWordsView();

  updateWordsScreenText(filteredCards.length);

  wordsDirectionButton.textContent =
    wordsDirection === "sv-pt" ?
    "Sueco → Português" :
    "Português → Sueco";

  wordsList.innerHTML = "";

  if (filteredCards.length === 0) {
    const emptyMessage = document.createElement("p");
    emptyMessage.className = "words-empty";
    emptyMessage.textContent = getEmptyWordsMessage();
    wordsList.appendChild(emptyMessage);
    return;
  }

  const sortedCards = sortCardsForWordsView(filteredCards);

  // Mantém "Revisar erros" como está hoje:
  // ordenado pela quantidade/taxa de erros.
  if (wordsViewMode === "wrong") {
    sortedCards.forEach((card) => {
      wordsList.appendChild(createWordItem(card));
    });

    return;
  }

  // Tela normal "Palavras":
  // agrupa alfabeticamente.
  const groupedCards = groupCardsByInitial(sortedCards);

  groupedCards.forEach((cardsInGroup, letter) => {
    wordsList.appendChild(
      createWordLetterGroup(letter, cardsInGroup)
    );
  });
}

function getCardsForWordsView() {
  const filteredCards = getFilteredCards();

  if (wordsViewMode !== "wrong") {
    return filteredCards;
  }

  return filteredCards.filter((card) => {
    const cardStats = getExistingCardStats(card.id);

    return cardStats && (cardStats.wrong || 0) > 0;
  });
}

function getExistingCardStats(cardId) {
  return stats[String(cardId)] || null;
}

function updateWordsScreenText(cardCount) {
  if (wordsViewMode === "wrong") {
    wordsScreenTitle.textContent = "Revisar erros";
    wordsScreenDescription.textContent =
      "Lista das palavras com erro salvo, respeitando os filtros escolhidos na configuração.";
    wordsCount.textContent =
      cardCount === 1 ?
      "1 palavra com erro" :
      `${cardCount} palavras com erro`;
    return;
  }

  wordsScreenTitle.textContent = "Palavras";
  wordsScreenDescription.textContent =
    "Lista gerada de acordo com os filtros escolhidos na configuração.";
  wordsCount.textContent =
    cardCount === 1 ?
    "1 palavra" :
    `${cardCount} palavras`;
}

function getEmptyWordsMessage() {
  if (wordsViewMode === "wrong") {
    return "Nenhuma palavra com erro encontrada com esses filtros.";
  }

  return "Nenhuma palavra encontrada com esses filtros.";
}

function getWordSortText(card) {
  if (wordsDirection === "sv-pt") {
    return card.term?.swedish || "";
  }

  return card.term?.portuguese || "";
}

function getWordInitial(card) {
  const text = getWordSortText(card).trim();

  if (!text) {
    return "#";
  }

  const locale =
    wordsDirection === "sv-pt" ?
    "sv-SE" :
    "pt-BR";

  let initial = Array.from(text)[0].toLocaleUpperCase(locale);

  // Em português, Á/Ã/Â ficam junto de A,
  // É/Ê junto de E etc.
  if (wordsDirection === "pt-sv") {
    initial = initial
      .normalize("NFD")
      .replace(/[\u0300-\u036f]/g, "");
  }

  return initial;
}

function groupCardsByInitial(cardList) {
  const groups = new Map();

  cardList.forEach((card) => {
    const letter = getWordInitial(card);

    if (!groups.has(letter)) {
      groups.set(letter, []);
    }

    groups.get(letter).push(card);
  });

  return groups;
}

function createWordLetterGroup(letter, cardList) {
  const details = document.createElement("details");
  details.className = "word-letter-group";
  details.dataset.letter = letter;

  // Reabre a letra depois que uma palavra individual
  // for expandida.
  details.open = expandedWordsLetter === letter;

  const summary = document.createElement("summary");
  summary.className = "word-letter-summary";

  const letterText = document.createElement("strong");
  letterText.textContent = letter;

  const count = document.createElement("span");
  count.className = "word-letter-count";
  count.textContent =
    `${cardList.length} ${cardList.length === 1 ? "palavra" : "palavras"}`;

  summary.append(letterText, count);

  const content = document.createElement("div");
  content.className = "word-letter-content";

  cardList.forEach((card) => {
    content.appendChild(createWordItem(card));
  });

  details.append(summary, content);

  details.addEventListener("toggle", () => {
    if (!details.open) {
      if (expandedWordsLetter === letter) {
        expandedWordsLetter = null;
      }

      return;
    }

    expandedWordsLetter = letter;

    // Fecha todas as outras letras.
    wordsList
      .querySelectorAll(".word-letter-group")
      .forEach((group) => {
        if (group !== details) {
          group.open = false;
        }
      });
  });

  return details;
}

function sortCardsForWordsView(cardList) {
  if (wordsViewMode === "wrong") {
    return [...cardList].sort((a, b) => {
      const aStats = getExistingCardStats(a.id) || {};
      const bStats = getExistingCardStats(b.id) || {};

      const wrongDifference = (bStats.wrong || 0) - (aStats.wrong || 0);

      if (wrongDifference !== 0) {
        return wrongDifference;
      }

      const aRate = getWrongRate(aStats);
      const bRate = getWrongRate(bStats);

      if (bRate !== aRate) {
        return bRate - aRate;
      }

      return getWordSortText(a).localeCompare(getWordSortText(b), "sv-SE");
    });
  }

  return [...cardList].sort((a, b) => {
    const aText = getWordSortText(a);
    const bText = getWordSortText(b);

    return aText.localeCompare(bText, "sv-SE");
  });
}

function getWrongRate(cardStats) {
  const correct = cardStats.correct || 0;
  const wrong = cardStats.wrong || 0;
  const total = correct + wrong;

  if (total === 0) {
    return 0;
  }

  return wrong / total;
}

function getWordStatsSummary(card) {
  const cardStats = getExistingCardStats(card.id) || {};
  const correct = cardStats.correct || 0;
  const wrong = cardStats.wrong || 0;
  const total = correct + wrong;
  const wrongRate = total === 0 ? 0 : Math.round((wrong / total) * 100);

  return `${wrong} erro${wrong === 1 ? "" : "s"} · ${correct} acerto${correct === 1 ? "" : "s"} · ${wrongRate}% de erro`;
}

function createWordItem(card) {
  const button = document.createElement("button");
  button.type = "button";
  button.className = "word-item";
  button.dataset.cardId = String(card.id);

  const main = document.createElement("div");
  main.className = "word-item-main";

  const primary = document.createElement("strong");
  primary.textContent = getWordPrimaryText(card);

  const secondary = document.createElement("span");
  secondary.textContent = getWordSecondaryText(card);

  main.append(primary, secondary);
  button.appendChild(main);

  if (String(card.id) === String(expandedWordCardId)) {
    button.appendChild(createWordDetails(card));
  }

  return button;
}

function createWordDetails(card) {
  const details = document.createElement("div");
  details.className = "word-details";

  const meta = document.createElement("p");
  meta.className = "word-meta";
  meta.textContent = getWordMetaText(card);
  details.appendChild(meta);

  const grammarDetails = createWordGrammarDetails(card);

  if (grammarDetails) {
    details.appendChild(grammarDetails);
  }

  const hasAudio = Boolean(card.media?.audio?.src);

  if (!hasAudio) {
    const noMedia = document.createElement("p");
    noMedia.className = "word-meta";
    noMedia.textContent = "Sem áudio para esta palavra.";
    details.appendChild(noMedia);
    return details;
  }

  const mediaActions = document.createElement("div");
  mediaActions.className = "word-media-actions";

  if (hasAudio) {
    const audioButton = document.createElement("button");
    audioButton.type = "button";
    audioButton.className = "word-audio-button";
    audioButton.dataset.audioSrc = card.media.audio.src;
    audioButton.textContent = "Ouvir áudio";

    mediaActions.appendChild(audioButton);
  }

  details.appendChild(mediaActions);

  return details;
}

function createWordGrammarDetails(card) {
  const grammar = card.grammar;

  if (!grammar) {
    return null;
  }

  const rows = [];

  if (grammar.type === "substantivo") {
    if (grammar.definiteSingular) {
      rows.push({
        label: "Singular definido",
        value: grammar.definiteSingular
      });
    }

    if (grammar.plural) {
      rows.push({
        label: "Plural indefinido",
        value: grammar.plural
      });
    }

    if (grammar.definitePlural) {
      rows.push({
        label: "Plural definido",
        value: grammar.definitePlural
      });
    }
  }

  if (grammar.type === "verbo") {
    if (grammar.infinitive) {
      rows.push({
        label: "Infinitivo",
        value: grammar.infinitive
      });
    }

    if (grammar.present) {
      rows.push({
        label: "Presente",
        value: grammar.present
      });
    }

    if (grammar.past) {
      rows.push({
        label: "Passado",
        value: grammar.past
      });
    }

    if (grammar.supine) {
      rows.push({
        label: "Supino",
        value: grammar.supine
      });
    }
  }

  if (rows.length === 0) {
    return null;
  }

  const container = document.createElement("div");
  container.className = "word-grammar-details";

  rows.forEach((row) => {
    const line = document.createElement("p");
    line.className = "word-grammar-row";

    const label = document.createElement("strong");
    label.textContent = `${row.label}:`;

    const value = document.createElement("span");
    value.textContent = row.value;

    line.append(label, value);
    container.appendChild(line);
  });

  return container;
}

function getWordSwedishListText(card) {
  const swedish = card.term?.swedish || "";

  if (!isNounWithGender(card)) {
    return swedish;
  }

  return `${swedish} · ${card.grammar.gender}`;
}

function getWordPrimaryText(card) {
  if (wordsDirection === "sv-pt") {
    return getWordSwedishListText(card);
  }

  return card.term.portuguese;
}

function getWordSecondaryText(card) {
  const translation = wordsDirection === "sv-pt" ?
    card.term.portuguese :
    getWordSwedishListText(card);

  if (wordsViewMode !== "wrong") {
    return translation;
  }

  return `${translation} · ${getWordStatsSummary(card)}`;
}

function getWordMetaText(card) {
  const parts = [];

  if (wordsViewMode === "wrong") {
    parts.push(getWordStatsSummary(card));
  }

  if (card.grammar?.type) {
    parts.push(`Tipo: ${card.grammar.type}`);
  }

  const themes = getCardThemes(card);

  if (themes.length > 0) {
    parts.push(`Tema: ${themes.join(", ")}`);
  }

  const source = getCardSource(card);

  if (source) {
    parts.push(`Origem: ${source}`);
  }

  if (source === "livro" && card.classification?.chapter) {
    parts.push(`Capítulo: ${card.classification.chapter}`);
  }

  if (source === "música" && getCardSourceTitle(card)) {
    parts.push(`Música: ${getCardSourceTitle(card)}`);
  }

  return parts.join(" · ");
}

function handleWordsListClick(event) {
  const audioButton = event.target.closest(".word-audio-button");

  if (audioButton) {
    event.stopPropagation();

    const audioSrc = audioButton.dataset.audioSrc;

    if (!audioSrc) {
      return;
    }

    const audio = new Audio(audioSrc);
    audio.play().catch((error) => {
      console.error("Erro ao tocar áudio da palavra:", error);
    });

    return;
  }

  const wordItem = event.target.closest(".word-item");

  if (!wordItem) {
    return;
  }

  const cardId = wordItem.dataset.cardId;

  if (String(expandedWordCardId) === String(cardId)) {
    expandedWordCardId = null;
  } else {
    expandedWordCardId = cardId;
  }

  renderWordsList();
}

function openExerciseScreen() {
  exerciseImportMessage.textContent = "";
  showScreen(exerciseScreen);
}

function backFromExerciseScreen() {
  showScreen(setupScreen);
}

function resetExerciseScreen() {
  currentExercise = null;
  exerciseFinished = false;
  exerciseErrorsCopyText = "";
  exerciseConsultedVocabulary = new Map();

  exerciseSourceInput.value = "";
  exerciseImportMessage.textContent = "";

  exerciseTitle.textContent = "";

  exerciseBlocks.replaceChildren();
  exerciseResultSummary.replaceChildren();

  exerciseImportPanel.classList.remove("hidden");
  exerciseRenderedPanel.classList.add("hidden");

  exerciseResultSummary.classList.add("hidden");

  finishExerciseButton.classList.remove("hidden");

  copyExerciseErrorsButton.classList.add("hidden");
  copyExerciseErrorsButton.textContent = "Copiar meus erros";

  newExerciseButton.classList.add("hidden");

  exerciseSourceInput.focus();
}

function generateExercise() {
  exerciseImportMessage.textContent = "";

  try {
    const parsedExercise = parseExerciseText(
      exerciseSourceInput.value
    );

    currentExercise = parsedExercise;
    exerciseFinished = false;
    exerciseErrorsCopyText = "";

    exerciseConsultedVocabulary = new Map();

    copyExerciseErrorsButton.classList.add("hidden");
    copyExerciseErrorsButton.textContent = "Copiar meus erros";

    renderExercise();

    exerciseImportPanel.classList.add("hidden");
    exerciseRenderedPanel.classList.remove("hidden");

    exerciseResultSummary.classList.add("hidden");

    finishExerciseButton.classList.remove("hidden");
    newExerciseButton.classList.add("hidden");

    window.scrollTo({
      top: 0,
      behavior: "smooth"
    });
  } catch (error) {
    console.error(error);

    exerciseImportMessage.textContent =
      error.message ||
      "Não foi possível interpretar o exercício.";
  }
}

function parseExerciseText(rawText) {
  const text = String(rawText || "")
    .replace(/\r\n?/g, "\n")
    .trim();

  if (!text) {
    throw new Error(
      "Cole um exercício antes de gerar."
    );
  }

  const lines = text.split("\n");

  if (lines[0].trim() !== "[EXERCICIO]") {
    throw new Error(
      'O exercício deve começar exatamente com "[EXERCICIO]".'
    );
  }

  let lastNonEmptyIndex = lines.length - 1;

  while (
    lastNonEmptyIndex >= 0 &&
    !lines[lastNonEmptyIndex].trim()
  ) {
    lastNonEmptyIndex--;
  }

  if (
    lastNonEmptyIndex < 0 ||
    lines[lastNonEmptyIndex].trim() !== "[FIM]"
  ) {
    throw new Error(
      'O exercício deve terminar exatamente com "[FIM]".'
    );
  }

  let title = "";
  const blocks = [];

  let questionNumber = 0;
  let i = 1;

  while (i <= lastNonEmptyIndex) {
    const line = lines[i].trim();

    if (!line) {
      i++;
      continue;
    }

    if (line === "[FIM]") {
      break;
    }

    if (line.startsWith("TITULO:")) {
      if (title) {
        throw new Error(
          `Linha ${i + 1}: existe mais de um campo TITULO.`
        );
      }

      title = line
        .slice("TITULO:".length)
        .trim();

      if (!title) {
        throw new Error(
          `Linha ${i + 1}: o campo TITULO está vazio.`
        );
      }

      i++;
      continue;
    }

    if (line === "[TEXTO]") {
      const startLine = i + 1;
      const contentLines = [];

      i++;

      while (
        i <= lastNonEmptyIndex &&
        !isExerciseStructuralMarker(
          lines[i].trim()
        )
      ) {
        contentLines.push(lines[i]);
        i++;
      }

      const content = contentLines
        .join("\n")
        .trim();

      if (!content) {
        throw new Error(
          `Linha ${startLine}: o bloco [TEXTO] está vazio.`
        );
      }

      blocks.push({
        kind: "text",
        text: content
      });

      continue;
    }

    if (line === "[QUESTAO]") {
      questionNumber++;

      const startLine = i + 1;
      const questionLines = [];

      i++;

      while (
        i <= lastNonEmptyIndex &&
        !isExerciseStructuralMarker(
          lines[i].trim()
        )
      ) {
        questionLines.push(lines[i]);
        i++;
      }

      blocks.push({
        kind: "question",

        question: parseExerciseQuestion(
          questionLines,
          questionNumber,
          startLine
        )
      });

      continue;
    }

    if (line === "[FALA]") {
      const startLine = i + 1;
      const speechLines = [];

      i++;

      while (
        i <= lastNonEmptyIndex &&
        !isExerciseStructuralMarker(
          lines[i].trim()
        )
      ) {
        speechLines.push(lines[i]);
        i++;
      }

      blocks.push({
        kind: "speech",

        items: parseExerciseSpeechBlock(
          speechLines,
          startLine
        )
      });

      continue;
    }

    throw new Error(
      `Linha ${i + 1}: conteúdo fora de um bloco reconhecido: "${line}".`
    );
  }

  if (!title) {
    throw new Error(
      "O exercício precisa de um campo TITULO:."
    );
  }

  const questions = blocks.filter(
    (block) => block.kind === "question"
  );

  if (questions.length === 0) {
    throw new Error(
      "O exercício precisa ter pelo menos uma [QUESTAO]."
    );
  }

  return {
    title,
    blocks,
    questionCount: questions.length
  };
}

function isExerciseStructuralMarker(line) {
  return (
    line === "[TEXTO]" ||
    line === "[QUESTAO]" ||
    line === "[FALA]" ||
    line === "[FIM]" ||
    line === "[EXERCICIO]"
  );
}

function parseExerciseSpeechBlock(
  lines,
  startLine
) {
  const meaningfulLines = lines
    .map((line, index) => ({
      text: line.trim(),
      lineNumber: startLine + index + 1
    }))
    .filter((item) => item.text);

  if (meaningfulLines.length === 0) {
    throw new Error(
      `Linha ${startLine}: o bloco [FALA] está vazio.`
    );
  }

  const items = [];
  let context = [];
  let i = 0;

  while (i < meaningfulLines.length) {
    const current = meaningfulLines[i];
    const line = current.text;

    if (line.startsWith("RESPOSTA:")) {
      throw new Error(
        `Linha ${current.lineNumber}: RESPOSTA: apareceu sem um VOCÊ: correspondente.`
      );
    }

    if (line.startsWith("VOCÊ:")) {
      const instruction = line
        .slice("VOCÊ:".length)
        .trim();

      if (!instruction) {
        throw new Error(
          `Linha ${current.lineNumber}: o campo VOCÊ: está vazio.`
        );
      }

      const answerLine =
        meaningfulLines[i + 1];

      if (
        !answerLine ||
        !answerLine.text.startsWith("RESPOSTA:")
      ) {
        throw new Error(
          `Linha ${current.lineNumber}: cada VOCÊ: deve ser seguido por RESPOSTA:.`
        );
      }

      const rawAnswer = answerLine.text
        .slice("RESPOSTA:".length)
        .trim();

      if (!rawAnswer) {
        throw new Error(
          `Linha ${answerLine.lineNumber}: o campo RESPOSTA está vazio.`
        );
      }

      const answers = rawAnswer
        .split("|")
        .map((answer) => answer.trim());

      if (
        answers.some(
          (answer) => !answer
        )
      ) {
        throw new Error(
          `Linha ${answerLine.lineNumber}: há uma resposta vazia ao redor do caractere |.`
        );
      }

      items.push({
        context,
        instruction,
        answers
      });

      context = [];
      i += 2;
      continue;
    }

    const dialogue =
      parseExerciseDialogueLine(line);

    if (!dialogue) {
      throw new Error(
        `Linha ${current.lineNumber}: dentro de [FALA], use "Nome: fala", "VOCÊ: instrução" ou "RESPOSTA: ...".`
      );
    }

    context.push(dialogue);
    i++;
  }

  if (context.length > 0) {
    throw new Error(
      "O bloco [FALA] terminou com uma fala do interlocutor sem um VOCÊ: correspondente."
    );
  }

  if (items.length === 0) {
    throw new Error(
      "O bloco [FALA] precisa ter pelo menos um VOCÊ: seguido de RESPOSTA:."
    );
  }

  return items;
}

function parseExerciseQuestion(
  lines,
  questionNumber,
  startLine
) {
  const meaningfulLines = lines
    .map((line, index) => ({
      text: line.trim(),
      lineNumber: startLine + index + 1
    }))
    .filter((item) => item.text);

  if (meaningfulLines.length === 0) {
    throw new Error(
      `Questão ${questionNumber}: a questão está vazia.`
    );
  }

  const firstLine = meaningfulLines[0];

  if (!firstLine.text.startsWith("TIPO:")) {
    throw new Error(
      `Questão ${questionNumber}: a primeira linha deve ser TIPO: MULTIPLA, TIPO: VF ou TIPO: ESCRITA.`
    );
  }

  const type = firstLine.text
    .slice("TIPO:".length)
    .trim()
    .toUpperCase();

  if (
    ![
      "MULTIPLA",
      "VF",
      "ESCRITA"
    ].includes(type)
  ) {
    throw new Error(
      `Questão ${questionNumber}: tipo "${type}" desconhecido. Use MULTIPLA, VF ou ESCRITA.`
    );
  }

  const promptLines = [];
  const options = [];

  let rawAnswer = "";
  let explanation = "";

  let answerFound = false;
  let explanationFound = false;

  for (
    const item of meaningfulLines.slice(1)
  ) {
    const line = item.text;

    if (line.startsWith("RESPOSTA:")) {
      if (answerFound) {
        throw new Error(
          `Questão ${questionNumber}: existe mais de um campo RESPOSTA.`
        );
      }

      rawAnswer = line
        .slice("RESPOSTA:".length)
        .trim();

      answerFound = true;

      if (!rawAnswer) {
        throw new Error(
          `Questão ${questionNumber}: o campo RESPOSTA está vazio.`
        );
      }

      continue;
    }

    if (line.startsWith("EXPLICACAO:")) {
      if (explanationFound) {
        throw new Error(
          `Questão ${questionNumber}: existe mais de um campo EXPLICACAO.`
        );
      }

      explanation = line
        .slice("EXPLICACAO:".length)
        .trim();

      explanationFound = true;

      if (!explanation) {
        throw new Error(
          `Questão ${questionNumber}: EXPLICACAO está vazia; remova a linha ou escreva a explicação.`
        );
      }

      continue;
    }

    if (type === "MULTIPLA") {
      const optionMatch =
        line.match(/^([A-Z])\)\s+(.+)$/);

      if (optionMatch) {
        options.push({
          letter: optionMatch[1],
          text: optionMatch[2].trim()
        });

        continue;
      }
    }

    if (
      answerFound ||
      explanationFound
    ) {
      throw new Error(
        `Questão ${questionNumber}: há conteúdo inesperado depois de RESPOSTA/EXPLICACAO: "${line}".`
      );
    }

    promptLines.push(line);
  }

  if (!answerFound) {
    throw new Error(
      `Questão ${questionNumber}: falta o campo RESPOSTA:.`
    );
  }

  if (type === "MULTIPLA") {
    const prompt = promptLines
      .join("\n")
      .trim();

    if (!prompt) {
      throw new Error(
        `Questão ${questionNumber}: falta o enunciado.`
      );
    }

    validateMultipleChoiceQuestion(
      options,
      rawAnswer,
      questionNumber
    );

    return {
      type,
      prompt,
      options,
      answer: rawAnswer.toUpperCase(),
      explanation
    };
  }

  if (type === "VF") {
    const prompt = promptLines
      .join("\n")
      .trim();

    if (!prompt) {
      throw new Error(
        `Questão ${questionNumber}: falta o enunciado.`
      );
    }

    const answer =
      rawAnswer.toUpperCase();

    if (
      answer !== "V" &&
      answer !== "F"
    ) {
      throw new Error(
        `Questão ${questionNumber}: questões VF aceitam apenas RESPOSTA: V ou RESPOSTA: F.`
      );
    }

    return {
      type,
      prompt,
      answer,
      explanation
    };
  }

  const groupedAnswers =
    parseExerciseGroupedWrittenAnswers(
      rawAnswer,
      questionNumber
    );

  if (groupedAnswers) {
    const groupedQuestion =
      parseExerciseGroupedWrittenPrompt(
        promptLines,
        groupedAnswers,
        questionNumber
      );

    return {
      type,
      prompt: groupedQuestion.prompt,
      subitems: groupedQuestion.subitems,
      explanation
    };
  }

  const prompt = promptLines
    .join("\n")
    .trim();

  if (!prompt) {
    throw new Error(
      `Questão ${questionNumber}: falta o enunciado.`
    );
  }

  const answerParts = rawAnswer
    .split("|")
    .map(
      (answer) => answer.trim()
    );

  if (
    answerParts.some(
      (answer) => !answer
    )
  ) {
    throw new Error(
      `Questão ${questionNumber}: há uma resposta vazia ao redor do caractere |.`
    );
  }

  return {
    type,
    prompt,
    answers: answerParts,
    explanation
  };
}

function parseExerciseGroupedWrittenAnswers(
  rawAnswer,
  questionNumber
) {
  if (!/^a\s*=/.test(rawAnswer)) {
    return null;
  }

  const markerRegex =
    /(?:^|;\s*)([a-z])\s*=/g;

  const matches = [
    ...rawAnswer.matchAll(markerRegex)
  ];

  if (matches.length < 2) {
    throw new Error(
      `Questão ${questionNumber}: uma questão escrita agrupada precisa ter pelo menos os subitens a) e b).`
    );
  }

  const entries = matches.map(
    (match, index) => {
      const letter = match[1];

      const valueStart =
        match.index + match[0].length;

      const valueEnd =
        index + 1 < matches.length ?
        matches[index + 1].index :
        rawAnswer.length;

      const rawValue = rawAnswer
        .slice(valueStart, valueEnd)
        .trim();

      const answers = rawValue
        .split("|")
        .map((answer) => answer.trim());

      if (
        answers.some(
          (answer) => !answer
        )
      ) {
        throw new Error(
          `Questão ${questionNumber}, subitem ${letter}): existe uma resposta vazia.`
        );
      }

      return {
        letter,
        answers
      };
    }
  );

  entries.forEach(
    (entry, index) => {
      const expectedLetter =
        String.fromCharCode(97 + index);

      if (
        entry.letter !== expectedLetter
      ) {
        throw new Error(
          `Questão ${questionNumber}: as respostas dos subitens devem ser sequenciais começando por a=. Esperado ${expectedLetter}=, encontrado ${entry.letter}=.`
        );
      }
    }
  );

  return entries;
}

function parseExerciseGroupedWrittenPrompt(
  promptLines,
  groupedAnswers,
  questionNumber
) {
  const instructionLines = [];
  const subitems = [];

  let currentSubitem = null;

  promptLines.forEach((line) => {
    const match =
      line.match(/^([a-z])\)\s*(.*)$/);

    if (match) {
      if (currentSubitem) {
        subitems.push(
          currentSubitem
        );
      }

      currentSubitem = {
        letter: match[1],
        lines: match[2] ? [match[2]] : []
      };

      return;
    }

    if (currentSubitem) {
      currentSubitem.lines.push(
        line
      );
    } else {
      instructionLines.push(
        line
      );
    }
  });

  if (currentSubitem) {
    subitems.push(
      currentSubitem
    );
  }

  const prompt = instructionLines
    .join("\n")
    .trim();

  if (!prompt) {
    throw new Error(
      `Questão ${questionNumber}: uma questão escrita agrupada precisa ter uma instrução geral antes do subitem a).`
    );
  }

  if (subitems.length < 2) {
    throw new Error(
      `Questão ${questionNumber}: RESPOSTA usa o formato agrupado, mas não foram encontrados pelo menos os subitens a) e b) no enunciado.`
    );
  }

  if (
    subitems.length !==
    groupedAnswers.length
  ) {
    throw new Error(
      `Questão ${questionNumber}: a quantidade de subitens (${subitems.length}) não corresponde à quantidade de respostas (${groupedAnswers.length}).`
    );
  }

  const parsedSubitems =
    subitems.map(
      (subitem, index) => {
        const expectedLetter =
          String.fromCharCode(
            97 + index
          );

        if (
          subitem.letter !==
          expectedLetter
        ) {
          throw new Error(
            `Questão ${questionNumber}: os subitens devem ser sequenciais começando por a). Esperado ${expectedLetter}), encontrado ${subitem.letter}).`
          );
        }

        const answerEntry =
          groupedAnswers[index];

        if (
          answerEntry.letter !==
          subitem.letter
        ) {
          throw new Error(
            `Questão ${questionNumber}: o subitem ${subitem.letter}) não corresponde à resposta ${answerEntry.letter}=.`
          );
        }

        const subitemPrompt =
          subitem.lines
          .join("\n")
          .trim();

        if (!subitemPrompt) {
          throw new Error(
            `Questão ${questionNumber}: o subitem ${subitem.letter}) está vazio.`
          );
        }

        return {
          letter: subitem.letter,
          prompt: subitemPrompt,
          answers: answerEntry.answers
        };
      }
    );

  return {
    prompt,
    subitems: parsedSubitems
  };
}

function validateMultipleChoiceQuestion(
  options,
  rawAnswer,
  questionNumber
) {
  if (options.length < 2) {
    throw new Error(
      `Questão ${questionNumber}: uma questão MULTIPLA precisa ter pelo menos duas alternativas.`
    );
  }

  options.forEach(
    (option, index) => {
      const expectedLetter =
        String.fromCharCode(
          65 + index
        );

      if (
        option.letter !==
        expectedLetter
      ) {
        throw new Error(
          `Questão ${questionNumber}: as alternativas devem ser sequenciais começando por A). Esperado ${expectedLetter}), encontrado ${option.letter}).`
        );
      }
    }
  );

  const answer =
    rawAnswer.toUpperCase();

  if (
    !/^[A-Z]$/.test(answer)
  ) {
    throw new Error(
      `Questão ${questionNumber}: RESPOSTA de múltipla escolha deve conter somente a letra da alternativa.`
    );
  }

  if (
    !options.some(
      (option) =>
      option.letter === answer
    )
  ) {
    throw new Error(
      `Questão ${questionNumber}: a resposta ${answer} não corresponde a nenhuma alternativa existente.`
    );
  }
}

function renderExercise() {
  exerciseTitle.textContent =
    currentExercise.title;

  exerciseBlocks.replaceChildren();
  exerciseResultSummary.replaceChildren();

  let textNumber = 0;
  let questionIndex = 0;

  currentExercise.blocks.forEach(
    (block) => {
      if (block.kind === "text") {
        textNumber++;

        exerciseBlocks.appendChild(
          createExerciseTextBlock(
            block.text,
            textNumber
          )
        );

        return;
      }

      if (block.kind === "speech") {
        exerciseBlocks.appendChild(
          createExerciseSpeechBlock(
            block.items
          )
        );

        return;
      }

      exerciseBlocks.appendChild(
        createExerciseQuestionBlock(
          block.question,
          questionIndex
        )
      );

      questionIndex++;
    }
  );
}

function createExerciseTextBlock(
  text,
  textNumber
) {
  const article =
    document.createElement("article");

  article.className =
    "exercise-reading-block";

  const label =
    document.createElement("p");

  label.className =
    "exercise-block-label";

  label.textContent =
    textNumber === 1 ?
    "Texto" :
    `Texto ${textNumber}`;

  const content =
    document.createElement("div");

  content.className =
    "exercise-reading-content";

  renderExerciseReadingContent(
    content,
    text
  );

  article.append(
    label,
    content
  );

  return article;
}

function renderExerciseReadingContent(
  container,
  text
) {
  container.replaceChildren();

  const lines = String(text || "")
    .replace(/\r\n?/g, "\n")
    .split("\n");

  let paragraphLines = [];

  const flushParagraph = () => {
    const paragraphText =
      paragraphLines
      .map((line) => line.trim())
      .filter(Boolean)
      .join(" ")
      .trim();

    paragraphLines = [];

    if (!paragraphText) {
      return;
    }

    const paragraph =
      document.createElement("p");

    paragraph.className =
      "exercise-reading-paragraph";

    renderExerciseInteractiveText(
      paragraph,
      paragraphText
    );

    container.appendChild(paragraph);
  };

  lines.forEach((rawLine) => {
    const line = rawLine.trim();

    if (!line) {
      flushParagraph();
      return;
    }

    const dialogue =
      parseExerciseDialogueLine(line);

    if (dialogue) {
      flushParagraph();

      container.appendChild(
        createExerciseDialogueLine(
          dialogue.speaker,
          dialogue.text
        )
      );

      return;
    }

    paragraphLines.push(line);
  });

  flushParagraph();
}

function parseExerciseDialogueLine(line) {
  const match = String(line || "").match(
    /^(\p{Lu}[\p{L}\p{M}'’.-]*(?:\s+\p{Lu}[\p{L}\p{M}'’.-]*){0,2}):\s*(.+)$/u
  );

  if (!match) {
    return null;
  }

  return {
    speaker: match[1].trim(),
    text: match[2].trim()
  };
}

function createExerciseDialogueLine(
  speakerName,
  dialogueText
) {
  const row =
    document.createElement("div");

  row.className =
    "exercise-dialogue-line";

  const speaker =
    document.createElement("strong");

  speaker.className =
    "exercise-dialogue-speaker";

  speaker.textContent =
    speakerName;

  const utterance =
    document.createElement("p");

  utterance.className =
    "exercise-dialogue-utterance";

  renderExerciseInteractiveText(
    utterance,
    dialogueText
  );

  row.append(
    speaker,
    utterance
  );

  return row;
}

function createExerciseSpeechBlock(items) {
  const article =
    document.createElement("article");

  article.className =
    "exercise-speech-block";

  const heading =
    document.createElement("div");

  heading.className =
    "exercise-speech-heading";

  const label =
    document.createElement("p");

  label.className =
    "exercise-block-label";

  label.textContent =
    "Prática oral";

  const intro =
    document.createElement("p");

  intro.className =
    "exercise-speech-intro";

  intro.textContent =
    "Fale em sueco e toque no campo azul para conferir respostas possíveis. Esta seção não vale nota.";

  heading.append(
    label,
    intro
  );

  const list =
    document.createElement("div");

  list.className =
    "exercise-speech-list";

  items.forEach((item, index) => {
    list.appendChild(
      createExerciseSpeechItem(
        item,
        index
      )
    );
  });

  article.append(
    heading,
    list
  );

  return article;
}

function createExerciseSpeechItem(
  item,
  itemIndex
) {
  const wrapper =
    document.createElement("div");

  wrapper.className =
    "exercise-speech-item";

  item.context.forEach(
    (dialogue) => {
      wrapper.appendChild(
        createExerciseDialogueLine(
          dialogue.speaker,
          dialogue.text
        )
      );
    }
  );

  const userRow =
    document.createElement("div");

  userRow.className =
    "exercise-speech-user-line";

  const speaker =
    document.createElement("strong");

  speaker.className =
    "exercise-dialogue-speaker";

  speaker.textContent = "Você";

  const responseArea =
    document.createElement("div");

  responseArea.className =
    "exercise-speech-response-area";

  const button =
    document.createElement("button");

  button.type = "button";

  button.className =
    "exercise-speech-prompt";

  button.textContent =
    item.instruction;

  const answerBox =
    document.createElement("div");

  answerBox.id =
    `exercise-speech-answer-${itemIndex}`;

  answerBox.className =
    "exercise-speech-answer hidden";

  const answerLabel =
    document.createElement("p");

  answerLabel.className =
    "exercise-speech-answer-label";

  answerLabel.textContent =
    item.answers.length === 1 ?
    "Resposta possível" :
    "Respostas possíveis";

  const answerList =
    document.createElement("div");

  answerList.className =
    "exercise-speech-answer-list";

  item.answers.forEach((answer) => {
    const paragraph =
      document.createElement("p");

    paragraph.className =
      "exercise-speech-answer-option";

    renderExerciseInteractiveText(
      paragraph,
      answer
    );

    answerList.appendChild(paragraph);
  });

  answerBox.append(
    answerLabel,
    answerList
  );

  button.setAttribute(
    "aria-expanded",
    "false"
  );

  button.setAttribute(
    "aria-controls",
    answerBox.id
  );

  button.setAttribute(
    "title",
    "Mostrar ou ocultar respostas possíveis"
  );

  button.addEventListener(
    "click",
    () => {
      const shouldShow =
        answerBox.classList.contains(
          "hidden"
        );

      answerBox.classList.toggle(
        "hidden",
        !shouldShow
      );

      button.classList.toggle(
        "revealed",
        shouldShow
      );

      button.setAttribute(
        "aria-expanded",
        String(shouldShow)
      );
    }
  );

  responseArea.append(
    button,
    answerBox
  );

  userRow.append(
    speaker,
    responseArea
  );

  wrapper.appendChild(userRow);

  return wrapper;
}

function buildExerciseVocabularyIndex() {
  exerciseVocabularyIndex = new Map();

  allCards.forEach((card) => {
    const swedish = String(
      card.term?.swedish || ""
    ).trim();

    if (!swedish) {
      return;
    }

    const words =
      extractExerciseSwedishWords(
        swedish
      );

    const isSingleWord =
      words.length === 1 &&
      normalizeExerciseVocabularyForm(
        words[0]
      ) ===
      normalizeExerciseVocabularyForm(
        swedish
      );

    /*
     * Mantém o registro original
     * exatamente como antes.
     */
    registerExerciseVocabularyForm(
      swedish, {
        source: "direct",
        priority: 30,
        card,
        formType: "base"
      }
    );

    /*
     * Cria somente uma chave alternativa
     * para busca, ignorando pontuação e
     * símbolos nas extremidades.
     *
     * Ex.:
     * "Säg till om jag kan göra något."
     * ->
     * "Säg till om jag kan göra något"
     */
    const searchableSwedish =
      swedish
      .replace(
        /^[\p{P}\p{S}\s]+|[\p{P}\p{S}\s]+$/gu,
        ""
      )
      .trim();

    /*
     * Só registra a forma alternativa
     * quando ela realmente for diferente.
     *
     * Prioridade 29:
     * - termo direto original = 30
     * - alias sem pontuação = 29
     * - flexões = 20
     * - palavra de expressão = 10
     */
    if (
      searchableSwedish &&
      normalizeExerciseVocabularyForm(
        searchableSwedish
      ) !==
      normalizeExerciseVocabularyForm(
        swedish
      )
    ) {
      registerExerciseVocabularyForm(
        searchableSwedish, {
          source: "direct",
          priority: 29,
          card,
          formType: "base"
        }
      );
    }

    /*
     * Comportamento antigo das expressões.
     */
    if (!isSingleWord) {
      words.forEach((word) => {
        registerExerciseVocabularyForm(
          word, {
            source: "expression",
            priority: 10,
            card,
            expression: swedish
          }
        );
      });
    }

    registerExerciseStoredInflections(
      card
    );

    /*
     * Mantém exatamente a lógica antiga
     * dos adjetivos.
     */
    registerExerciseRegularAdjectiveForms(
      card,
      isSingleWord ? swedish : ""
    );
  });
}

function registerExerciseVocabularyForm(
  form,
  entry
) {
  const key =
    normalizeExerciseVocabularyForm(form);

  if (!key) {
    return;
  }

  const current =
    exerciseVocabularyIndex.get(key);

  if (
    !current ||
    entry.priority > current.priority
  ) {
    exerciseVocabularyIndex.set(
      key,
      entry
    );
  }
}

function registerExerciseStoredInflections(
  card
) {
  const grammar = card.grammar || {};

  const forms = [];

  if (grammar.type === "substantivo") {
    forms.push(
      ["plural", "plural"],
      [
        "definiteSingular",
        "definiteSingular"
      ],
      [
        "definitePlural",
        "definitePlural"
      ]
    );
  }

  if (grammar.type === "verbo") {
    forms.push(
      ["infinitive", "infinitive"],
      ["present", "present"],
      ["past", "past"],
      ["supine", "supine"]
    );
  }

  forms.forEach(([field, formType]) => {
    const form = grammar[field];

    if (!form) {
      return;
    }

    registerExerciseVocabularyForm(
      form, {
        source: "inflection",
        priority: 20,
        card,
        formType
      }
    );
  });
}

function registerExerciseRegularAdjectiveForms(
  card,
  swedish
) {
  if (
    card.grammar?.type !== "adjetivo" ||
    !swedish
  ) {
    return;
  }

  const base =
    normalizeExerciseVocabularyForm(
      swedish
    );

  if (
    !base ||
    !/^[\p{L}\p{M}]+$/u.test(base)
  ) {
    return;
  }

  // Reconhece apenas o padrão regular mais seguro:
  // base + t para ett
  // base + a para plural/definido.
  //
  // Formas com alteração de radical ou outras
  // irregularidades devem continuar explícitas
  // no JSON para reconhecimento perfeito.

  if (
    /[aeiouyåäö]$/u.test(base) ||
    /(?:al|el|er|en|d|t)$/u.test(base)
  ) {
    return;
  }

  registerExerciseVocabularyForm(
    `${base}t`, {
      source: "inflection",
      priority: 20,
      card,
      formType: "adjectiveNeuter"
    }
  );

  registerExerciseVocabularyForm(
    `${base}a`, {
      source: "inflection",
      priority: 20,
      card,
      formType: "adjectivePlural"
    }
  );
}

function extractExerciseSwedishWords(text) {
  return String(text || "").match(
    /[\p{L}\p{M}]+(?:['’\-][\p{L}\p{M}]+)*/gu
  ) || [];
}

function normalizeExerciseVocabularyForm(
  value
) {
  return String(value || "")
    .normalize("NFC")
    .toLocaleLowerCase("sv-SE")
    .trim();
}

function renderExerciseInteractiveText(
  container,
  text
) {
  container.replaceChildren();

  const source = String(text || "");

  const wordRegex =
    /[\p{L}\p{M}]+(?:['’\-][\p{L}\p{M}]+)*/gu;

  const words = Array.from(
    source.matchAll(wordRegex),
    (match) => ({
      word: match[0],
      start: match.index || 0,
      end: (match.index || 0) +
        match[0].length
    })
  );

  if (words.length === 0) {
    container.appendChild(
      document.createTextNode(source)
    );

    return;
  }

  /*
   * Descobre o maior número de palavras
   * existente em uma entrada do índice.
   *
   * Exemplos:
   *
   * "bok" = 1
   * "i måndags" = 2
   * "lade sig" = 2
   */
  let maxVocabularyWords = 1;

  for (
    const key of
      exerciseVocabularyIndex.keys()
  ) {
    const count =
      extractExerciseSwedishWords(
        key
      ).length;

    if (count > maxVocabularyWords) {
      maxVocabularyWords = count;
    }
  }

  let cursor = 0;
  let wordIndex = 0;

  while (wordIndex < words.length) {
    const current =
      words[wordIndex];

    /*
     * Preserva pontuação, espaços etc.
     * existentes antes da palavra atual.
     */
    if (current.start > cursor) {
      container.appendChild(
        document.createTextNode(
          source.slice(
            cursor,
            current.start
          )
        )
      );
    }

    let matchedEntry = null;
    let matchedEndIndex =
      wordIndex;

    const maximumEndIndex =
      Math.min(
        words.length - 1,
        wordIndex +
        maxVocabularyWords -
        1
      );

    /*
     * Procura primeiro a expressão
     * mais longa possível.
     *
     * Assim:
     *
     * "i måndags"
     *
     * é procurado antes de:
     *
     * "i"
     *
     * e:
     *
     * "lade sig"
     *
     * é procurado antes de:
     *
     * "lade"
     */
    for (
      let endIndex =
        maximumEndIndex; endIndex >= wordIndex; endIndex--
    ) {
      let validSequence = true;

      /*
       * Entre as palavras da expressão
       * pode existir apenas espaço.
       *
       * "i måndags"     -> válido
       * "i, måndags"    -> inválido
       * "lade sig"      -> válido
       * "lade. sig"     -> inválido
       */
      for (
        let separatorIndex =
          wordIndex; separatorIndex < endIndex; separatorIndex++
      ) {
        const separator =
          source.slice(
            words[separatorIndex].end,
            words[
              separatorIndex + 1
            ].start
          );

        const onlyWhitespace =
          /^\s+$/u.test(separator);

        const containsLineBreak =
          /[\r\n]/u.test(separator);

        if (
          !onlyWhitespace ||
          containsLineBreak
        ) {
          validSequence = false;
          break;
        }
      }

      if (!validSequence) {
        continue;
      }

      const candidate =
        words
        .slice(
          wordIndex,
          endIndex + 1
        )
        .map((item) => item.word)
        .join(" ");

      const entry =
        exerciseVocabularyIndex.get(
          normalizeExerciseVocabularyForm(
            candidate
          )
        );

      if (!entry) {
        continue;
      }

      matchedEntry = entry;
      matchedEndIndex =
        endIndex;

      /*
       * Estamos indo da expressão
       * maior para a menor.
       * A primeira encontrada vence.
       */
      break;
    }

    if (matchedEntry) {
      const finalWord =
        words[matchedEndIndex];

      /*
       * Usa o texto original para manter
       * maiúsculas/minúsculas exatamente
       * como aparecem no exercício.
       */
      const visibleText =
        source.slice(
          current.start,
          finalWord.end
        );

      container.appendChild(
        createExerciseVocabularyToken(
          visibleText,
          matchedEntry
        )
      );

      cursor = finalWord.end;

      wordIndex =
        matchedEndIndex + 1;

      continue;
    }

    /*
     * Nada correspondente no JSON.
     */
    container.appendChild(
      document.createTextNode(
        current.word
      )
    );

    cursor = current.end;
    wordIndex++;
  }

  /*
   * Preserva o que vier depois
   * da última palavra.
   */
  if (cursor < source.length) {
    container.appendChild(
      document.createTextNode(
        source.slice(cursor)
      )
    );
  }
}

function createExerciseVocabularyToken(
  word,
  entry
) {
  const wrapper =
    document.createElement("span");

  wrapper.className =
    "exercise-vocab-token";

  const button =
    document.createElement("button");

  button.type = "button";

  button.className =
    "exercise-vocab-word";

  button.textContent = word;

  button.title =
    "Consultar no words.json";

  button.setAttribute(
    "aria-expanded",
    "false"
  );

  const details =
    createExerciseVocabularyDetails(
      word,
      entry
    );

  details.classList.add("hidden");

  button.addEventListener(
    "click",
    () => {
      const shouldOpen =
        details.classList.contains(
          "hidden"
        );

      closeExpandedExerciseVocabulary();

      if (!shouldOpen) {
        return;
      }

      registerExerciseVocabularyConsultation(
        word,
        entry
      );

      details.classList.remove(
        "hidden"
      );

      button.setAttribute(
        "aria-expanded",
        "true"
      );

      expandedExerciseVocabulary = {
        button,
        details
      };
    }
  );

  wrapper.append(
    button,
    details
  );

  return wrapper;
}

function registerExerciseVocabularyConsultation(
  word,
  entry
) {
  const visibleForm =
    String(word || "").trim();

  if (
    !visibleForm ||
    !entry ||
    !entry.card
  ) {
    return;
  }

  const card = entry.card;

  const baseSwedish =
    String(
      card.term?.swedish ||
      visibleForm
    ).trim();

  const portuguese =
    String(
      card.term?.portuguese || ""
    ).trim();

  /*
   * Palavras extraídas de uma expressão
   * precisam permanecer independentes.
   *
   * Exemplo:
   *
   * "gärna" encontrada dentro de
   * "Jag kommer gärna"
   *
   * não deve ser agrupada como se a
   * expressão inteira tivesse sido clicada.
   */
  const isExpressionWord =
    entry.source === "expression";

  const cardIdentifier =
    card.id != null ?
    String(card.id) :
    normalizeExerciseVocabularyForm(
      baseSwedish
    );

  const key =
    isExpressionWord ? [
      "expression",
      cardIdentifier,
      normalizeExerciseVocabularyForm(
        visibleForm
      )
    ].join(":") : [
      "card",
      cardIdentifier
    ].join(":");

  const current =
    exerciseConsultedVocabulary.get(
      key
    );

  if (current) {
    current.count++;

    current.forms.add(
      visibleForm
    );

    return;
  }

  exerciseConsultedVocabulary.set(
    key, {
      swedish: isExpressionWord ?
        visibleForm : baseSwedish,

      portuguese,

      count: 1,

      forms: new Set([
        visibleForm
      ]),

      expression: isExpressionWord ?
        String(
          entry.expression ||
          baseSwedish
        ).trim() : ""
    }
  );
}

function closeExpandedExerciseVocabulary() {
  if (!expandedExerciseVocabulary) {
    return;
  }

  expandedExerciseVocabulary.details
    .classList.add("hidden");

  expandedExerciseVocabulary.button
    .setAttribute(
      "aria-expanded",
      "false"
    );

  expandedExerciseVocabulary = null;
}

function createExerciseVocabularyDetails(
  word,
  entry
) {
  const details =
    document.createElement("span");

  details.className =
    "exercise-vocab-details";

  details.setAttribute(
    "role",
    "note"
  );

  const title =
    document.createElement("span");

  title.className =
    "exercise-vocab-detail-title";

  title.textContent = word;

  details.appendChild(title);

  if (entry.source === "expression") {
    const info =
      document.createElement("span");

    info.className =
      "exercise-vocab-detail-meta";

    info.textContent =
      "Sem tradução individual cadastrada. Palavra encontrada em uma expressão do JSON:";

    const expression =
      document.createElement("span");

    expression.className =
      "exercise-vocab-detail-expression";

    expression.textContent =
      entry.expression;

    const translation =
      document.createElement("span");

    translation.className =
      "exercise-vocab-detail-translation";

    translation.textContent =
      entry.card.term?.portuguese || "";

    details.append(
      info,
      expression,
      translation
    );

    return details;
  }

  const translation =
    getExercisePortugueseForEntry(
      entry
    );

  if (translation) {
    const translationLine =
      document.createElement("span");

    translationLine.className =
      "exercise-vocab-detail-translation";

    translationLine.textContent =
      translation;

    details.appendChild(
      translationLine
    );
  }

  const formDescription =
    getExerciseVocabularyFormDescription(
      entry
    );

  if (formDescription) {
    const formLine =
      document.createElement("span");

    formLine.className =
      "exercise-vocab-detail-meta";

    formLine.textContent =
      formDescription;

    details.appendChild(formLine);
  }

  return details;
}

function getExercisePortugueseForEntry(
  entry
) {
  const term =
    entry.card?.term || {};

  switch (entry.formType) {
    case "infinitive":
      return (
        term.portugueseInfinitive ||
        term.portuguese ||
        ""
      );

    case "present":
      return (
        term.portuguesePresent ||
        term.portuguese ||
        ""
      );

    case "past":
      return (
        term.portuguesePast ||
        term.portuguese ||
        ""
      );

    case "supine":
      return (
        term.portugueseSupine ||
        term.portuguese ||
        ""
      );

    case "plural":
    case "definitePlural":
      return (
        term.portuguesePlural ||
        term.portuguese ||
        ""
      );

    default:
      return (
        term.portuguese ||
        ""
      );
  }
}

function getExerciseVocabularyFormDescription(
  entry
) {
  if (
    entry.source !== "inflection"
  ) {
    return "";
  }

  const base = String(
    entry.card?.term?.swedish || ""
  ).trim();

  const labels = {
    infinitive: "Infinitivo",

    present: "Presente",

    past: "Pretérito",

    supine: "Supino",

    plural: "Plural",

    definiteSingular: "Singular definido",

    definitePlural: "Plural definido",

    adjectiveNeuter: "Forma do adjetivo com substantivo ett",

    adjectivePlural: "Forma do adjetivo no plural/definido"
  };

  const label =
    labels[entry.formType];

  if (!label || !base) {
    return "";
  }

  return `${label} de ${base}`;
}

function createExerciseQuestionBlock(
  question,
  questionIndex
) {
  const article =
    document.createElement("article");

  article.className =
    "exercise-question-card";

  article.dataset.exerciseQuestionIndex =
    String(questionIndex);

  const number =
    document.createElement("p");

  number.className =
    "exercise-question-number";

  number.textContent =
    `Questão ${questionIndex + 1}`;

  const prompt =
    document.createElement("p");

  prompt.className =
    "exercise-question-prompt";

  prompt.textContent =
    question.prompt;

  article.append(
    number,
    prompt
  );

  if (
    question.type === "MULTIPLA"
  ) {
    article.appendChild(
      createExerciseMultipleOptions(
        question,
        questionIndex
      )
    );
  } else if (
    question.type === "VF"
  ) {
    article.appendChild(
      createExerciseTrueFalseOptions(
        questionIndex
      )
    );
  } else if (
    Array.isArray(question.subitems) &&
    question.subitems.length > 0
  ) {
    article.appendChild(
      createExerciseGroupedWrittenInputs(
        question,
        questionIndex
      )
    );
  } else {
    article.appendChild(
      createExerciseWrittenInput(
        questionIndex
      )
    );
  }

  const feedback =
    document.createElement("div");

  feedback.className =
    "exercise-feedback hidden";

  feedback.dataset.exerciseFeedback =
    String(questionIndex);

  article.appendChild(feedback);

  return article;
}

function createExerciseMultipleOptions(
  question,
  questionIndex
) {
  const wrapper =
    document.createElement("div");

  wrapper.className =
    "exercise-options";

  question.options.forEach(
    (option) => {
      wrapper.appendChild(
        createExerciseOption(
          `exercise-question-${questionIndex}`,
          option.letter,
          `${option.letter}) ${option.text}`
        )
      );
    }
  );

  return wrapper;
}

function createExerciseTrueFalseOptions(
  questionIndex
) {
  const wrapper =
    document.createElement("div");

  wrapper.className =
    "exercise-options exercise-options-vf";

  wrapper.append(
    createExerciseOption(
      `exercise-question-${questionIndex}`,
      "V",
      "Verdadeiro"
    ),

    createExerciseOption(
      `exercise-question-${questionIndex}`,
      "F",
      "Falso"
    )
  );

  return wrapper;
}

function createExerciseOption(
  name,
  value,
  labelText
) {
  const label =
    document.createElement("label");

  label.className =
    "exercise-option";

  const input =
    document.createElement("input");

  input.type = "radio";
  input.name = name;
  input.value = value;

  const text =
    document.createElement("span");

  text.textContent = labelText;

  label.append(
    input,
    text
  );

  return label;
}

function createExerciseGroupedWrittenInputs(
  question,
  questionIndex
) {
  const wrapper =
    document.createElement("div");

  wrapper.className =
    "exercise-written-group";

  question.subitems.forEach(
    (subitem) => {
      const item =
        document.createElement("div");

      item.className =
        "exercise-written-subitem";

      const label =
        document.createElement("label");

      label.className =
        "exercise-written-subitem-prompt";

      const inputId =
        `exercise-written-${questionIndex}-${subitem.letter}`;

      label.htmlFor = inputId;

      const letter =
        document.createElement("strong");

      letter.className =
        "exercise-written-subitem-letter";

      letter.textContent =
        `${subitem.letter})`;

      const text =
        document.createElement("span");

      text.className =
        "exercise-written-subitem-text";

      text.textContent =
        subitem.prompt;

      label.append(
        letter,
        text
      );

      const input =
        createExerciseWrittenInput(
          questionIndex,
          subitem.letter
        );

      input.id = inputId;

      item.append(
        label,
        input
      );

      wrapper.appendChild(item);
    }
  );

  return wrapper;
}

function createExerciseWrittenInput(
  questionIndex,
  subitemLetter = null
) {
  const input =
    document.createElement("input");

  input.className =
    "exercise-written-input";

  input.type = "text";

  // Evita interferência do teclado/navegador
  input.setAttribute("autocomplete", "off");
  input.setAttribute("autocorrect", "off");
  input.setAttribute("autocapitalize", "none");
  input.setAttribute("spellcheck", "false");
  input.setAttribute("writingsuggestions", "false");

  // Indica que o conteúdo esperado é sueco
  input.setAttribute("lang", "sv");

  input.setAttribute("inputmode", "text");
  input.setAttribute("enterkeyhint", "done");

  input.dataset.exerciseWrittenInput =
    String(questionIndex);

  if (subitemLetter) {
    input.dataset.exerciseSubitem =
      subitemLetter;

    input.placeholder =
      `Resposta ${subitemLetter})...`;
  } else {
    input.placeholder =
      "Digite sua resposta...";
  }

  return input;
}

function finishExercise() {
  if (
    !currentExercise ||
    exerciseFinished
  ) {
    return;
  }

  const questions =
    currentExercise.blocks
    .filter(
      (block) =>
      block.kind === "question"
    )
    .map(
      (block) =>
      block.question
    );

  let totalScore = 0;

  const counts = {
    correct: 0,
    partial: 0,
    wrong: 0,
    unanswered: 0
  };

  const exerciseErrors = [];

  questions.forEach(
    (question, questionIndex) => {
      const article =
        exerciseBlocks.querySelector(
          `[data-exercise-question-index="${questionIndex}"]`
        );

      const userAnswer =
        getExerciseUserAnswer(
          question,
          questionIndex,
          article
        );

      const result =
        gradeExerciseQuestion(
          question,
          userAnswer
        );

      totalScore += result.score;

      counts[result.status]++;

      if (result.status !== "correct") {
        exerciseErrors.push({
          question,
          questionIndex,
          userAnswer,
          result
        });
      }

      renderExerciseQuestionFeedback(
        article,
        question,
        userAnswer,
        result
      );

      article
        .querySelectorAll("input")
        .forEach(
          (input) => {
            input.disabled = true;
          }
        );
    }
  );

  exerciseErrorsCopyText =
    buildExerciseErrorsCopyText(
      exerciseErrors
    );

  exerciseFinished = true;

  finishExerciseButton.classList.add(
    "hidden"
  );

  copyExerciseErrorsButton.textContent =
    "Copiar meus erros";

  copyExerciseErrorsButton.classList.toggle(
    "hidden",
    !exerciseErrorsCopyText
  );

  newExerciseButton.classList.remove(
    "hidden"
  );

  renderExerciseResultSummary(
    totalScore,
    questions.length,
    counts
  );

  exerciseResultSummary.classList.remove(
    "hidden"
  );

  exerciseResultSummary.scrollIntoView({
    behavior: "smooth",
    block: "start"
  });
}

function getExerciseReviewCorrection(
  userAnswer,
  expectedAnswer
) {
  const userWords =
    String(userAnswer || "")
    .trim()
    .split(/\s+/)
    .map((word) => normalizeAnswer(word))
    .filter(Boolean);

  const expectedWordEntries =
    String(expectedAnswer || "")
    .trim()
    .split(/\s+/)
    .map((word) => ({
      original: word,
      normalized: normalizeAnswer(word)
    }))
    .filter((item) => item.normalized);

  const expectedWords =
    expectedWordEntries.map(
      (item) => item.normalized
    );

  if (expectedWords.length === 0) {
    return "";
  }

  if (userWords.length === 0) {
    return String(expectedAnswer || "").trim();
  }

  const comparison =
    getExerciseWordComparison(
      userWords,
      expectedWords
    );

  const correctedIndexes = [
    ...new Set(
      comparison.operations
      .filter(
        (operation) =>
        (
          operation.type === "substitute" ||
          operation.type === "insert"
        ) &&
        Number.isInteger(
          operation.expectedIndex
        )
      )
      .map(
        (operation) =>
        operation.expectedIndex
      )
    )
  ].sort((a, b) => a - b);

  if (correctedIndexes.length === 0) {
    return String(expectedAnswer || "").trim();
  }

  return correctedIndexes
    .map(
      (index) =>
      expectedWordEntries[index]?.original ||
      expectedWords[index]
    )
    .filter(Boolean)
    .join(" ");
}


function getExerciseConsultedWordsForCopy() {
  return [
    ...new Set(
      Array.from(
        exerciseConsultedVocabulary.values()
      )
      .flatMap((item) => {
        const forms =
          item.forms instanceof Set ?
          Array.from(item.forms) :
          [];

        if (forms.length > 0) {
          return forms;
        }

        return item.swedish ?
          [item.swedish] :
          [];
      })
      .map(
        (word) =>
        String(word || "").trim()
      )
      .filter(Boolean)
    )
  ].sort(
    (a, b) =>
    a.localeCompare(b, "sv")
  );
}


function buildExerciseErrorsCopyText(
  errorEntries
) {
  const errors =
    Array.isArray(errorEntries) ?
    errorEntries :
    [];

  const consultedWords =
    getExerciseConsultedWordsForCopy();

  if (
    errors.length === 0 &&
    consultedWords.length === 0
  ) {
    return "";
  }

  const lines = [
    "Revisão",
    ""
  ];

  errors.forEach(
    ({
      question,
      userAnswer,
      result
    }) => {
      /*
       * Questões escritas com subitens:
       * não copia número da questão,
       * letra, instrução, resposta do usuário
       * nem explicação do erro.
       */
      if (
        question.type === "ESCRITA" &&
        Array.isArray(question.subitems) &&
        question.subitems.length > 0
      ) {
        result.subitemResults
          .filter(
            (subitemResult) =>
            subitemResult.status !==
            "correct"
          )
          .forEach(
            (subitemResult) => {
              const subitem =
                question.subitems.find(
                  (item) =>
                  item.letter ===
                  subitemResult.letter
                );

              if (subitem?.prompt) {
                lines.push(
                  subitem.prompt
                );
              }

              lines.push(
                `-> ${getExerciseReviewCorrection(
                  subitemResult.userAnswer,
                  subitemResult.expectedAnswer
                )}`
              );

              lines.push("");
            }
          );

        return;
      }

      /*
       * Questões sem subitens.
       */
      if (question.prompt) {
        lines.push(
          question.prompt
        );
      }

      const correction =
        question.type === "ESCRITA" ?
        getExerciseReviewCorrection(
          userAnswer,
          result.expectedAnswer
        ) :
        result.expectedAnswer;

      lines.push(
        `-> ${correction}`
      );

      lines.push("");
    }
  );

  /*
   * Adiciona somente as palavras que
   * foram realmente consultadas.
   * Sem tradução e sem quantidade.
   */
  if (consultedWords.length > 0) {
    lines.push(
      "Palavras consultadas:"
    );

    consultedWords.forEach(
      (word) => {
        lines.push(word);
      }
    );
  }

  return lines
    .join("\n")
    .trim();
}

async function copyExerciseErrors() {
  if (!exerciseErrorsCopyText) {
    return;
  }

  try {
    await navigator.clipboard.writeText(
      exerciseErrorsCopyText
    );

    copyExerciseErrorsButton.textContent =
      "Erros copiados";
  } catch (error) {
    console.error(
      "Não foi possível copiar pelo Clipboard API:",
      error
    );

    const textarea =
      document.createElement("textarea");

    textarea.value =
      exerciseErrorsCopyText;

    textarea.setAttribute(
      "readonly",
      ""
    );

    textarea.style.position = "fixed";
    textarea.style.opacity = "0";

    document.body.appendChild(
      textarea
    );

    textarea.select();

    try {
      const copied =
        document.execCommand("copy");

      if (!copied) {
        throw new Error(
          "Comando de cópia não disponível."
        );
      }

      copyExerciseErrorsButton.textContent =
        "Erros copiados";
    } catch (fallbackError) {
      console.error(
        "Também não foi possível copiar pelo método alternativo:",
        fallbackError
      );

      copyExerciseErrorsButton.textContent =
        "Não foi possível copiar";
    } finally {
      textarea.remove();
    }
  }
}

function getExerciseUserAnswer(
  question,
  questionIndex,
  article
) {
  if (
    question.type === "ESCRITA" &&
    Array.isArray(question.subitems) &&
    question.subitems.length > 0
  ) {
    const answers = {};

    question.subitems.forEach(
      (subitem) => {
        const input =
          article.querySelector(
            `[data-exercise-written-input="${questionIndex}"][data-exercise-subitem="${subitem.letter}"]`
          );

        answers[subitem.letter] =
          input ?
          input.value.trim() :
          "";
      }
    );

    return answers;
  }

  if (
    question.type === "ESCRITA"
  ) {
    const input =
      article.querySelector(
        `[data-exercise-written-input="${questionIndex}"]`
      );

    return input ?
      input.value.trim() :
      "";
  }

  const checked =
    article.querySelector(
      `input[name="exercise-question-${questionIndex}"]:checked`
    );

  return checked ?
    checked.value :
    "";
}

function gradeExerciseQuestion(
  question,
  userAnswer
) {
  if (
    question.type === "ESCRITA" &&
    Array.isArray(question.subitems) &&
    question.subitems.length > 0
  ) {
    return gradeExerciseGroupedWrittenQuestion(
      question,
      userAnswer
    );
  }

  const isEmptyWrittenAnswer =
    question.type === "ESCRITA" &&
    !normalizeAnswer(userAnswer);

  if (
    !userAnswer ||
    isEmptyWrittenAnswer
  ) {
    return {
      score: 0,
      status: "unanswered",

      expectedAnswer: getExerciseExpectedAnswerText(
        question
      ),

      feedback: [
        "Questão não respondida."
      ]
    };
  }

  if (
    question.type === "MULTIPLA" ||
    question.type === "VF"
  ) {
    const isCorrect =
      userAnswer ===
      question.answer;

    return {
      score: isCorrect ? 1 : 0,

      status: isCorrect ?
        "correct" : "wrong",

      expectedAnswer: getExerciseExpectedAnswerText(
        question
      ),

      feedback: [
        isCorrect ?
        "Resposta correta." :
        "Resposta incorreta."
      ]
    };
  }

  return compareExerciseWrittenAnswer(
    userAnswer,
    question.answers
  );
}

function gradeExerciseGroupedWrittenQuestion(
  question,
  userAnswers
) {
  const answers =
    userAnswers &&
    typeof userAnswers === "object" ?
    userAnswers : {};

  const subitemResults =
    question.subitems.map(
      (subitem) => {
        const userAnswer =
          answers[subitem.letter] || "";

        if (!normalizeAnswer(userAnswer)) {
          return {
            letter: subitem.letter,
            userAnswer,
            score: 0,
            status: "unanswered",
            expectedAnswer: subitem.answers[0] || "",
            feedback: [
              "Subitem não respondido."
            ]
          };
        }

        const result =
          compareExerciseWrittenAnswer(
            userAnswer,
            subitem.answers
          );

        return {
          ...result,
          letter: subitem.letter,
          userAnswer
        };
      }
    );

  const score =
    subitemResults.reduce(
      (sum, result) =>
      sum + result.score,
      0
    ) / subitemResults.length;

  const allUnanswered =
    subitemResults.every(
      (result) =>
      result.status === "unanswered"
    );

  let status = "wrong";

  if (allUnanswered) {
    status = "unanswered";
  } else if (score === 1) {
    status = "correct";
  } else if (score > 0) {
    status = "partial";
  }

  return {
    score,
    status,
    expectedAnswer: question.subitems
      .map(
        (subitem) =>
        `${subitem.letter}) ${subitem.answers[0] || ""}`
      )
      .join("; "),
    feedback: [],
    subitemResults
  };
}

function getExerciseExpectedAnswerText(
  question
) {
  if (
    question.type === "MULTIPLA"
  ) {
    const option =
      question.options.find(
        (item) =>
        item.letter ===
        question.answer
      );

    return option ?
      `${option.letter}) ${option.text}` :
      question.answer;
  }

  if (
    question.type === "VF"
  ) {
    return question.answer === "V" ?
      "Verdadeiro" :
      "Falso";
  }

  return question.answers[0] || "";
}

function getExerciseDisplayedUserAnswer(
  question,
  userAnswer
) {
  if (!userAnswer) {
    return "Não respondida";
  }

  if (
    question.type === "MULTIPLA"
  ) {
    const option =
      question.options.find(
        (item) =>
        item.letter === userAnswer
      );

    return option ?
      `${option.letter}) ${option.text}` :
      userAnswer;
  }

  if (
    question.type === "VF"
  ) {
    return userAnswer === "V" ?
      "Verdadeiro" :
      "Falso";
  }

  return userAnswer;
}

function compareExerciseWrittenAnswer(
  userAnswer,
  acceptedAnswers
) {
  const normalizedUserAnswer =
    normalizeAnswer(userAnswer);

  const candidates =
    acceptedAnswers.map(
      (acceptedAnswer) => {
        const normalizedExpectedAnswer =
          normalizeAnswer(
            acceptedAnswer
          );

        if (
          normalizedUserAnswer ===
          normalizedExpectedAnswer
        ) {
          return {
            score: 1,

            status: "correct",

            expectedAnswer: acceptedAnswer,

            normalizedExpectedAnswer,

            feedback: [
              "Resposta correta."
            ],

            charDistance: 0
          };
        }

        return evaluateExerciseWrittenCandidate(
          normalizedUserAnswer,
          normalizedExpectedAnswer,
          acceptedAnswer
        );
      }
    );

  candidates.sort(
    (a, b) => {
      const similarityA =
        a.similarityScore ?? a.score;

      const similarityB =
        b.similarityScore ?? b.score;

      if (
        similarityB !== similarityA
      ) {
        return (
          similarityB -
          similarityA
        );
      }

      return (
        a.charDistance -
        b.charDistance
      );
    }
  );

  return candidates[0];
}

function evaluateExerciseWrittenCandidate(
  normalizedUserAnswer,
  normalizedExpectedAnswer,
  originalExpectedAnswer
) {
  const charComparison =
    getDamerauLevenshteinComparison(
      normalizedUserAnswer,
      normalizedExpectedAnswer
    );

  const maxCharLength =
    Math.max(
      normalizedUserAnswer.length,
      normalizedExpectedAnswer.length,
      1
    );

  const orthographyScore =
    clamp(
      1 -
      charComparison.distance /
      maxCharLength,
      0,
      1
    );

  const userWords =
    normalizedUserAnswer
    .split(/\s+/)
    .filter(Boolean);

  const expectedWords =
    normalizedExpectedAnswer
    .split(/\s+/)
    .filter(Boolean);

  const wordComparison =
    getExerciseWordComparison(
      userWords,
      expectedWords
    );

  const maxWordLength =
    Math.max(
      userWords.length,
      expectedWords.length,
      1
    );

  const contentScore =
    clamp(
      1 -
      wordComparison.distance /
      maxWordLength,
      0,
      1
    );

  // Usado SOMENTE para diagnosticar
  // quão próxima a resposta ficou.
  const similarityScore =
    clamp(
      contentScore * 0.7 +
      orthographyScore * 0.3,
      0,
      0.9999
    );

  return {
    score: 0,

    similarityScore,

    status: "wrong",

    expectedAnswer: originalExpectedAnswer,

    normalizedExpectedAnswer,

    charDistance: charComparison.distance,

    feedback: buildExerciseWrittenFeedback(
      wordComparison.operations
    )
  };
}

function getExerciseWordComparison(
  sourceWords,
  targetWords
) {
  const sourceLength =
    sourceWords.length;

  const targetLength =
    targetWords.length;

  const distances =
    Array.from({
        length: sourceLength + 1
      },
      () =>
      Array(
        targetLength + 1
      ).fill(0)
    );

  const steps =
    Array.from({
        length: sourceLength + 1
      },
      () =>
      Array(
        targetLength + 1
      ).fill(null)
    );

  for (
    let i = 1; i <= sourceLength; i++
  ) {
    distances[i][0] = i;

    steps[i][0] = {
      type: "delete"
    };
  }

  for (
    let j = 1; j <= targetLength; j++
  ) {
    distances[0][j] = j;

    steps[0][j] = {
      type: "insert"
    };
  }

  for (
    let i = 1; i <= sourceLength; i++
  ) {
    for (
      let j = 1; j <= targetLength; j++
    ) {
      const userWord =
        sourceWords[i - 1];

      const expectedWord =
        targetWords[j - 1];

      const similarity =
        getExerciseWordSimilarity(
          userWord,
          expectedWord
        );

      let substitutionCost = 0;

      if (
        userWord !== expectedWord
      ) {
        substitutionCost =
          similarity >= 0.6 ?
          (1 - similarity) * 0.6 :
          1;
      }

      let bestDistance =
        distances[i - 1][j - 1] +
        substitutionCost;

      let bestStep = {
        type: substitutionCost === 0 ?
          "match" : "substitute",

        similarity
      };

      const deleteDistance =
        distances[i - 1][j] + 1;

      if (
        deleteDistance <
        bestDistance
      ) {
        bestDistance =
          deleteDistance;

        bestStep = {
          type: "delete"
        };
      }

      const insertDistance =
        distances[i][j - 1] + 1;

      if (
        insertDistance <
        bestDistance
      ) {
        bestDistance =
          insertDistance;

        bestStep = {
          type: "insert"
        };
      }

      distances[i][j] =
        bestDistance;

      steps[i][j] =
        bestStep;
    }
  }

  const operations = [];

  let i = sourceLength;
  let j = targetLength;

  while (
    i > 0 ||
    j > 0
  ) {
    const step =
      steps[i][j];

    if (!step) {
      break;
    }

    if (
      step.type === "match"
    ) {
      i--;
      j--;

      continue;
    }

    if (
      step.type === "substitute"
    ) {
      operations.push({
        type: "substitute",

        userWord: sourceWords[i - 1],

        expectedWord: targetWords[j - 1],

        userIndex: i - 1,

        expectedIndex: j - 1,

        similarity: step.similarity
      });

      i--;
      j--;

      continue;
    }

    if (
      step.type === "delete"
    ) {
      operations.push({
        type: "delete",

        userWord: sourceWords[i - 1],

        userIndex: i - 1,

        expectedIndex: j
      });

      i--;

      continue;
    }

    if (
      step.type === "insert"
    ) {
      operations.push({
        type: "insert",

        expectedWord: targetWords[j - 1],

        userIndex: i,

        expectedIndex: j - 1
      });

      j--;
    }
  }

  operations.reverse();

  return {
    distance: distances[
      sourceLength
    ][
      targetLength
    ],

    operations
  };
}

function getExerciseWordSimilarity(
  userWord,
  expectedWord
) {
  if (
    userWord === expectedWord
  ) {
    return 1;
  }

  const comparison =
    getDamerauLevenshteinComparison(
      userWord,
      expectedWord
    );

  const maxLength =
    Math.max(
      userWord.length,
      expectedWord.length,
      1
    );

  return clamp(
    1 -
    comparison.distance /
    maxLength,
    0,
    1
  );
}

function buildExerciseWrittenFeedback(
  operations
) {
  if (!operations.length) {
    return [
      "Resposta correta."
    ];
  }

  const missingWords =
    operations
    .filter(
      (operation) =>
      operation.type ===
      "insert"
    )
    .map(
      (operation) =>
      operation.expectedWord
    );

  const extraWords =
    operations
    .filter(
      (operation) =>
      operation.type ===
      "delete"
    )
    .map(
      (operation) =>
      operation.userWord
    );

  const feedback = [];

  if (
    missingWords.length === 1
  ) {
    feedback.push(
      `Faltou a palavra "${missingWords[0]}".`
    );
  } else if (
    missingWords.length > 1
  ) {
    feedback.push(
      `Faltaram as palavras ${formatExerciseQuotedList(missingWords)}.`
    );
  }

  if (
    extraWords.length === 1
  ) {
    feedback.push(
      `Há uma palavra extra: "${extraWords[0]}".`
    );
  } else if (
    extraWords.length > 1
  ) {
    feedback.push(
      `Há palavras extras: ${formatExerciseQuotedList(extraWords)}.`
    );
  }

  operations
    .filter(
      (operation) =>
      operation.type ===
      "substitute"
    )
    .forEach(
      (operation) => {
        if (
          operation.similarity >=
          0.6
        ) {
          feedback.push(
            ...buildExerciseWordSpellingFeedback(
              operation.userWord,
              operation.expectedWord,
              operation.expectedIndex +
              1
            )
          );
        } else {
          feedback.push(
            `Na ${getOrdinalWord(operation.expectedIndex + 1)} palavra, você escreveu "${operation.userWord}", mas o esperado era "${operation.expectedWord}".`
          );
        }
      }
    );

  if (
    feedback.length === 0
  ) {
    return [
      "A resposta ficou diferente da forma esperada."
    ];
  }

  return feedback.slice(0, 6);
}

function buildExerciseWordSpellingFeedback(
  userWord,
  expectedWord,
  wordNumber
) {
  const comparison =
    getDamerauLevenshteinComparison(
      userWord,
      expectedWord
    );

  const prefix =
    `Na ${getOrdinalWord(wordNumber)} palavra`;

  if (
    comparison.operations.length >
    3
  ) {
    return [
      `${prefix}, você escreveu "${userWord}", mas o correto é "${expectedWord}".`
    ];
  }

  return comparison.operations.map(
    (operation) => {
      if (
        operation.type ===
        "substitute"
      ) {
        const part =
          getWordPart(
            clamp(
              operation.position -
              1,
              0,
              expectedWord.length -
              1
            ),

            expectedWord.length
          );

        return `${prefix}, ${part}, você escreveu "${operation.userChar}", mas o correto é "${operation.expectedChar}".`;
      }

      if (
        operation.type ===
        "insert"
      ) {
        const part =
          getWordPart(
            clamp(
              operation.position -
              1,
              0,
              expectedWord.length -
              1
            ),

            expectedWord.length
          );

        return `${prefix}, faltou a letra "${operation.expectedChar}" ${part}.`;
      }

      if (
        operation.type ===
        "delete"
      ) {
        const part =
          getWordPart(
            clamp(
              operation.position -
              1,
              0,
              userWord.length -
              1
            ),

            userWord.length
          );

        return `${prefix}, há uma letra extra "${operation.userChar}" ${part}.`;
      }

      if (
        operation.type ===
        "transpose"
      ) {
        return `${prefix}, algumas letras parecem estar invertidas.`;
      }

      return `${prefix}, a grafia está diferente.`;
    }
  );
}

function formatExerciseQuotedList(
  items
) {
  const quoted =
    items.map(
      (item) => `"${item}"`
    );

  if (
    quoted.length <= 1
  ) {
    return quoted[0] || "";
  }

  if (
    quoted.length === 2
  ) {
    return (
      `${quoted[0]} e ${quoted[1]}`
    );
  }

  return (
    `${quoted
      .slice(0, -1)
      .join(", ")} e ${quoted[quoted.length - 1]}`
  );
}

function renderExerciseQuestionFeedback(
  article,
  question,
  userAnswer,
  result
) {
  const feedbackBox =
    article.querySelector(
      "[data-exercise-feedback]"
    );

  feedbackBox.replaceChildren();

  feedbackBox.className =
    `exercise-feedback exercise-feedback-${result.status}`;

  if (
    question.type === "ESCRITA" &&
    Array.isArray(question.subitems) &&
    question.subitems.length > 0
  ) {
    renderExerciseGroupedWrittenFeedback(
      feedbackBox,
      question,
      result
    );

    return;
  }

  const heading =
    document.createElement(
      "strong"
    );

  heading.className =
    "exercise-feedback-title";

  heading.textContent =
    getExerciseStatusLabel(
      result.status,
      result.score,
      question.type
    );

  const userLine =
    document.createElement("p");

  userLine.textContent =
    `Sua resposta: ${getExerciseDisplayedUserAnswer(question, userAnswer)}`;

  const expectedLine =
    document.createElement("p");

  expectedLine.textContent =
    `Resposta esperada: ${result.expectedAnswer}`;

  feedbackBox.append(
    heading,
    userLine,
    expectedLine
  );

  result.feedback.forEach(
    (message) => {
      if (
        message ===
        "Resposta correta."
      ) {
        return;
      }

      const paragraph =
        document.createElement(
          "p"
        );

      paragraph.className =
        "exercise-feedback-detail";

      paragraph.textContent =
        message;

      feedbackBox.appendChild(
        paragraph
      );
    }
  );

  appendExerciseExplanation(
    feedbackBox,
    question
  );
}

function renderExerciseGroupedWrittenFeedback(
  feedbackBox,
  question,
  result
) {
  const heading =
    document.createElement(
      "strong"
    );

  heading.className =
    "exercise-feedback-title";

  heading.textContent =
    getExerciseStatusLabel(
      result.status,
      result.score,
      "ESCRITA"
    );

  const list =
    document.createElement("div");

  list.className =
    "exercise-subitem-feedback-list";

  result.subitemResults.forEach(
    (subitemResult) => {
      const item =
        document.createElement("div");

      item.className =
        `exercise-subitem-feedback exercise-subitem-feedback-${subitemResult.status}`;

      const title =
        document.createElement(
          "strong"
        );

      title.className =
        "exercise-subitem-feedback-title";

      title.textContent =
        `${subitemResult.letter}) ${getExerciseStatusLabel(
          subitemResult.status,
          subitemResult.score,
          "ESCRITA"
        )}`;

      const userLine =
        document.createElement("p");

      userLine.textContent =
        `Sua resposta: ${subitemResult.userAnswer || "Não respondida"}`;

      const expectedLine =
        document.createElement("p");

      expectedLine.textContent =
        `Resposta esperada: ${subitemResult.expectedAnswer}`;

      item.append(
        title,
        userLine,
        expectedLine
      );

      subitemResult.feedback.forEach(
        (message) => {
          if (
            message === "Resposta correta." ||
            message === "Subitem não respondido."
          ) {
            return;
          }

          const detail =
            document.createElement("p");

          detail.className =
            "exercise-feedback-detail";

          detail.textContent =
            message;

          item.appendChild(detail);
        }
      );

      list.appendChild(item);
    }
  );

  feedbackBox.append(
    heading,
    list
  );

  appendExerciseExplanation(
    feedbackBox,
    question
  );
}

function appendExerciseExplanation(
  feedbackBox,
  question
) {
  if (!question.explanation) {
    return;
  }

  const explanation =
    document.createElement("p");

  explanation.className =
    "exercise-explanation";

  explanation.textContent =
    `Explicação: ${question.explanation}`;

  feedbackBox.appendChild(
    explanation
  );
}

function getExerciseStatusLabel(
  status,
  score,
  questionType
) {
  if (
    status === "correct"
  ) {
    return "✅ Correta";
  }

  if (
    status === "unanswered"
  ) {
    return "⚪ Não respondida — Nota 0,0";
  }

  if (
    status === "partial"
  ) {
    return (
      `🟡 Parcialmente correta — Nota ${formatExerciseGrade(score)}`
    );
  }

  return "Incorreto";
}

function createExerciseConsultedVocabularySummary() {
  const section =
    document.createElement("section");

  section.className =
    "exercise-consulted-vocabulary";

  const title =
    document.createElement("h4");

  title.textContent =
    "Palavras consultadas";

  section.appendChild(title);

  const items =
    Array.from(
      exerciseConsultedVocabulary.values()
    ).sort(
      (a, b) => {
        if (b.count !== a.count) {
          return b.count - a.count;
        }

        return a.swedish.localeCompare(
          b.swedish,
          "sv"
        );
      }
    );

  if (items.length === 0) {
    const empty =
      document.createElement("p");

    empty.className =
      "exercise-consulted-empty";

    empty.textContent =
      "Nenhuma palavra consultada.";

    section.appendChild(empty);

    return section;
  }

  const list =
    document.createElement("div");

  list.className =
    "exercise-consulted-list";

  items.forEach((item) => {
    const row =
      document.createElement("div");

    row.className =
      "exercise-consulted-item";

    const header =
      document.createElement("div");

    header.className =
      "exercise-consulted-item-header";

    const word =
      document.createElement("strong");

    word.textContent =
      item.swedish;

    const count =
      document.createElement("span");

    count.textContent =
      item.count === 1 ?
      "1 consulta" :
      `${item.count} consultas`;

    header.append(
      word,
      count
    );

    row.appendChild(header);

    if (item.portuguese) {
      const translation =
        document.createElement("p");

      translation.className =
        "exercise-consulted-translation";

      translation.textContent =
        item.portuguese;

      row.appendChild(
        translation
      );
    }

    const forms =
      Array.from(item.forms);

    const normalizedBase =
      normalizeExerciseVocabularyForm(
        item.swedish
      );

    const shouldShowForms =
      forms.length > 1 ||
      (
        forms.length === 1 &&
        normalizeExerciseVocabularyForm(
          forms[0]
        ) !== normalizedBase
      );

    if (
      shouldShowForms &&
      !item.expression
    ) {
      const formInfo =
        document.createElement("p");

      formInfo.className =
        "exercise-consulted-meta";

      formInfo.textContent =
        `Forma${forms.length > 1 ? "s" : ""} consultada${forms.length > 1 ? "s" : ""}: ${forms.join(", ")}`;

      row.appendChild(
        formInfo
      );
    }

    if (item.expression) {
      const expressionInfo =
        document.createElement("p");

      expressionInfo.className =
        "exercise-consulted-meta";

      expressionInfo.textContent =
        `Encontrada em: ${item.expression}`;

      row.appendChild(
        expressionInfo
      );
    }

    list.appendChild(row);
  });

  section.appendChild(list);

  return section;
}

function renderExerciseResultSummary(
  totalScore,
  questionCount,
  counts
) {
  exerciseResultSummary.replaceChildren();

  const title =
    document.createElement("h3");

  title.textContent =
    "Resultado do exercício";

  const finalScore =
    questionCount > 0 ?
    (
      totalScore /
      questionCount
    ) :
    0;

  const score =
    document.createElement("p");

  score.className =
    "exercise-result-score";

  score.textContent =
    `Nota: ${formatExerciseGrade(finalScore)}`;

  const grid =
    document.createElement("div");

  grid.className =
    "exercise-result-grid";

  appendExerciseResultMetric(
    grid,
    "Corretas",
    counts.correct
  );

  appendExerciseResultMetric(
    grid,
    "Parciais",
    counts.partial
  );

  appendExerciseResultMetric(
    grid,
    "Incorretas",
    counts.wrong
  );

  appendExerciseResultMetric(
    grid,
    "Não respondidas",
    counts.unanswered
  );

  const consultedVocabulary =
    createExerciseConsultedVocabularySummary();

  exerciseResultSummary.append(
    title,
    consultedVocabulary,
    score,
    grid
  );
}

function appendExerciseResultMetric(
  container,
  label,
  value
) {
  const item =
    document.createElement("div");

  const labelElement =
    document.createElement("span");

  labelElement.textContent =
    label;

  const valueElement =
    document.createElement(
      "strong"
    );

  valueElement.textContent =
    String(value);

  item.append(
    labelElement,
    valueElement
  );

  container.appendChild(
    item
  );
}

function formatExercisePoints(
  value
) {
  return value.toLocaleString(
    "pt-BR", {
      minimumFractionDigits: 0,
      maximumFractionDigits: 2
    }
  );
}

function formatExerciseGrade(
  score
) {
  const normalizedScore =
    clamp(score, 0, 1);

  if (normalizedScore === 1) {
    return "10";
  }

  const grade =
    Math.min(
      9.9,
      Math.round(
        normalizedScore * 100
      ) / 10
    );

  return grade.toLocaleString(
    "pt-BR", {
      minimumFractionDigits: 1,
      maximumFractionDigits: 1
    }
  );
}

function resetStats() {
  const filteredCards = getFilteredCards();

  if (filteredCards.length === 0) {
    setupMessage.textContent = "Nenhuma palavra encontrada com os filtros atuais.";
    return;
  }

  const shouldReset = window.confirm(
    `Tem certeza que deseja zerar o progresso de ${filteredCards.length} palavra(s) dos filtros atuais? As demais palavras permanecerão inalteradas.`
  );

  if (!shouldReset) {
    return;
  }

  filteredCards.forEach((card) => {
    delete stats[String(card.id)];
  });

  if (Object.keys(stats).length === 0) {
    localStorage.removeItem(STORAGE_KEY);
  } else {
    saveStats();
  }

  correctCount = 0;
  wrongCount = 0;
  sessionAnswers = [];

  setupMessage.textContent = `Progresso zerado para ${filteredCards.length} palavra(s) dos filtros atuais. As demais palavras foram mantidas.`;
}

const audioPreloadCache = new Map();

function preloadAudioForCard(card) {
  const audioSrc = card?.media?.audio?.src;

  if (!audioSrc || audioPreloadCache.has(audioSrc)) {
    return;
  }

  const audio = new Audio();
  audio.preload = "auto";
  audio.src = audioSrc;
  audio.load();

  audioPreloadCache.set(audioSrc, audio);
}

function preloadAudioForCards(cardList, limit = 6) {
  cardList
    .filter((card) => card?.media?.audio?.src)
    .slice(0, limit)
    .forEach(preloadAudioForCard);
}

flashcard.addEventListener("click", revealAnswer);

flashcard.addEventListener("keydown", (event) => {
  if (event.key === "Enter" || event.key === " ") {
    event.preventDefault();
    revealAnswer();
  }
});

checkAnswerButton.addEventListener("click", checkWrittenAnswer);
nextWriteCardButton.addEventListener("click", nextWrittenCard);

answerInput.addEventListener("keydown", (event) => {
  if (event.key === "Enter") {
    checkWrittenAnswer();
  }
});

playAudioButton.addEventListener("click", (event) => {
  event.stopPropagation();
  cardAudio.currentTime = 0;
  cardAudio.play();
});

showSwedishFromAudioButton.addEventListener("click", (event) => {
  event.stopPropagation();
  showSwedishFromAudio();
});

themeToggleButton.addEventListener("click", toggleTheme);

directionSelect.addEventListener("change", updateModeUI);
answerModeSelect.addEventListener("change", updateModeUI);

startSessionButton.addEventListener("click", startSession);
backToSetupButton.addEventListener("click", backToSetup);

pronunciationRulesButton.addEventListener("click", openPronunciationRules);
backFromPronunciationButton.addEventListener("click", backFromPronunciationRules);
grammarButton.addEventListener("click", openGrammarScreen);
backFromGrammarButton.addEventListener("click", backFromGrammarScreen);

pluralRulesButton.addEventListener("click", openPluralRules);
backFromPluralButton.addEventListener("click", backFromPluralRules);

verbFormsButton.addEventListener("click", openVerbForms);
backFromVerbFormsButton.addEventListener("click", backFromVerbForms);

pronunciationRulesList.addEventListener("click", playPronunciationExample);

resetStatsButton.addEventListener("click", resetStats);
wordsButton.addEventListener("click", () => openWordsScreen("all"));
wrongWordsButton.addEventListener("click", () => openWordsScreen("wrong"));
backFromWordsButton.addEventListener("click", backFromWordsScreen);
wordsDirectionButton.addEventListener("click", toggleWordsDirection);
wordsList.addEventListener("click", handleWordsListClick);

exerciseButton.addEventListener(
  "click",
  openExerciseScreen
);

backFromExerciseButton.addEventListener(
  "click",
  backFromExerciseScreen
);

generateExerciseButton.addEventListener(
  "click",
  generateExercise
);

finishExerciseButton.addEventListener(
  "click",
  finishExercise
);

copyExerciseErrorsButton.addEventListener(
  "click",
  copyExerciseErrors
);

newExerciseButton.addEventListener(
  "click",
  resetExerciseScreen
);

finishSessionButton.addEventListener("click", showSummary);
repeatSessionButton.addEventListener("click", repeatSession);
newSessionButton.addEventListener("click", backToSetup);

correctButton.addEventListener("click", () => markAnswer(true));
wrongButton.addEventListener("click", () => markAnswer(false));

newWordsToggleButton.addEventListener("click", toggleNewWordsMode);

holdCurrentWordButton.addEventListener(
  "click",
  holdCurrentNewWord
);

restoreHeldWordsButton.addEventListener(
  "click",
  restoreHeldNewWords
);

sourceFilterGroup.addEventListener("change", () => {
  setupMessage.textContent = "";
  updateSourceSpecificFilters();
});

pronunciationButton.addEventListener(
  "click",
  playCurrentCardPronunciation
);

applySavedTheme();

loadHeldNewWords();

updateNewWordsModeUI();
updateHeldNewWordsUI();

loadCards();